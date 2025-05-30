import os
import numpy as np
from mayavi import mlab
from tvtk.tools import visual
from scipy.interpolate import Rbf
from scipy.interpolate import LinearNDInterpolator

os.environ['ETS_TOOLKIT'] = 'qt'
os.environ['QT_API'] = 'pyqt5'

from mayavi.core.ui.api import MayaviScene, SceneEditor, \
                MlabSceneModel
from traits.api import HasTraits, Range, Instance, \
        on_trait_change

def Arrow_From_A_to_B(x1, y1, z1, x2, y2, z2):
    """
    Create an arrow from A to B given by coordinates x1, y1, z1 for A and x2, y2, z2 for B
    """
    ar1=visual.arrow(x=x1, y=y1, z=z1)
    ar1.color = (0.5, 0.5, 0.5)

    arrow_length=np.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
    ar1.length_cone=0.02*30/arrow_length
    ar1.radius_shaft = 0.015*30/arrow_length
    ar1.radius_cone = 0.04*30/arrow_length
    ar1.actor.scale=[arrow_length, arrow_length, arrow_length]
    ar1.pos = ar1.pos/arrow_length
    ar1.axis = [x2-x1, y2-y1, z2-z1]
    return ar1

def create_folder_unless_it_exists(path):
    """
    Create folder if it does not exist. If it exists, do nothing.

    Parameters
    ----------
    path: str
        Path to the folder to be created or checked.
    """
    if not os.path.exists(path):
        os.makedirs(path)


def plot_3d_dataset_as_cube(x_raw, y_raw, z_raw, k_raw, substance_titles = ('Alcohol', 'HBr', 'Temperature'),
                            npoints=30, sparse_npoints=4, rbf_epsilon=0.01, rbf_smooth=0.001,
                            interpolator_choice='rbf', data_for_spheres='raw', colormap='blue-red', contours=5,
                            colorbar_title="Yield", rbf_function="multiquadric", axes_ticks_format='%.2f', axes_font_factor=0.83,
                            contour_opacity=1, single_size_of_points=False,
                            forced_kmax=None, dont_mlabshow=False, reorient_callable=None, savefig=None,
                            transparent=False):

    xs0 = (x_raw - np.min(x_raw)) / (np.max(x_raw) - np.min(x_raw))
    ys0 = (y_raw - np.min(y_raw)) / (np.max(y_raw) - np.min(y_raw))
    zs0 = (z_raw - np.min(z_raw)) / (np.max(z_raw) - np.min(z_raw))
    ks0 = k_raw
    max_xs0 = np.max(xs0)
    max_ys0 = np.max(ys0)
    max_zs0 = np.max(zs0)
    max_ks0 = np.max(ks0)

    if forced_kmax is not None:
        max_ks0 = forced_kmax

    npoints = npoints * 1j
    xnew, ynew, znew = np.mgrid[np.min(xs0):np.max(xs0):npoints,
                       np.min(ys0):np.max(ys0):npoints,
                       np.min(zs0):np.max(zs0):npoints]

    # construct sparse points for spheres
    sparse_npoints_complex = 3*[sparse_npoints * 1j]
    xnew_sparse, ynew_sparse, znew_sparse = np.mgrid[np.min(xs0):np.max(xs0):sparse_npoints_complex[0],
                                            np.min(ys0):np.max(ys0):sparse_npoints_complex[1],
                                            np.min(zs0):np.max(zs0):sparse_npoints_complex[2]]

    if interpolator_choice == 'rbf':
        rbf4 = Rbf(xs0, ys0, zs0, ks0, epsilon=rbf_epsilon, smooth=rbf_smooth, function=rbf_function)
        wnew = rbf4(xnew, ynew, znew)
        if data_for_spheres == 'interpolated':
            sparse_ks = rbf4(xnew_sparse, ynew_sparse, znew_sparse)
    elif interpolator_choice == 'linear':
        interp_here = LinearNDInterpolator((xs0, ys0, zs0), ks0)
        wnew = interp_here((xnew, ynew, znew))
        if data_for_spheres == 'interpolated':
            sparse_ks = interp_here((xnew_sparse, ynew_sparse, znew_sparse))
    else:
        raise ValueError('interpolator_choice must be either "rbf" or "linear"')

    if data_for_spheres == 'raw':
        xs, ys, zs, ks = xs0, ys0, zs0, ks0
    elif data_for_spheres == 'interpolated':
        xs, ys, zs, ks = xnew_sparse, ynew_sparse, znew_sparse, sparse_ks
    else:
        raise ValueError('data_for_spheres must be either "raw" or "interpolated"')

    mlab.figure(size=(1024, 1224), bgcolor=(1, 1, 1), fgcolor=(0.2, 0.2, 0.2))
    plot = mlab.contour3d(xnew, ynew, znew, wnew, extent=[np.min(xs0), max_xs0,
                                                          np.min(ys0), max_ys0,
                                                          np.min(zs0), max_zs0],
                          contours=contours, opacity=contour_opacity, vmin=0, vmax=max_ks0, colormap=colormap,
                          transparent=transparent)
    plot.actor.actor.property.ambient = 0.0
    ax1 = mlab.axes(color=(0.5, 0.5, 0.5), nb_labels=sparse_npoints, ranges=[np.min(x_raw), np.max(x_raw),
                                                                np.min(y_raw), np.max(y_raw),
                                                                np.min(z_raw), np.max(z_raw)])
    mlab.xlabel(f'{substance_titles[0]}')
    mlab.ylabel(f'{substance_titles[1]}')
    mlab.zlabel(f'{substance_titles[2]}')
    mlab.outline(plot)
    cb = mlab.colorbar(object=plot, title=colorbar_title, orientation='horizontal', nb_labels=5)
    cb.scalar_bar.unconstrained_font_size = True
    cb.label_text_property.font_size = 19
    ax1.axes.font_factor = axes_font_factor
    ax1.axes.label_format= axes_ticks_format
    ax1.axes.corner_offset = 0.05
    if forced_kmax is not None:
        ks_scaled = (ks - np.min(ks)) / (np.max(ks) - np.min(ks))
        plot_points = mlab.points3d(xs, ys, zs, ks_scaled, vmin=0, vmax=np.max((forced_kmax - np.min(ks)) / (np.max(ks) - np.min(ks))),
                                    resolution=16, scale_factor=0.1, colormap=colormap)
    else:
        ks_scaled = (ks - np.min(ks)) / (np.max(ks) - np.min(ks))
        plot_points = mlab.points3d(xs, ys, zs, ks_scaled, vmin=0, vmax=np.max((ks0 - np.min(ks)) / (np.max(ks) - np.min(ks))),
                                    resolution=16, scale_factor=0.1, colormap=colormap)
    # plot_points = mlab.points3d(xs, ys, zs, ks_scaled, vmin=0, vmax=np.max(ks_scaled),
    #                             resolution=16, scale_factor=0.1, colormap=colormap)

    # get current mayavi scene
    if reorient_callable is not None:
        scene = mlab.get_engine().scenes[0]
        reorient_callable(scene)

    if savefig is not None:
        mlab.savefig(savefig)

    if not dont_mlabshow:
        mlab.show()


if __name__ == '__main__':
    # Initial generation of the data
    npoints = 7j
    x_raw, y_raw, z_raw = np.mgrid[1:5:npoints,
                          10:50:npoints,
                          100:500:npoints]
    # flatten all arrays
    x_raw = x_raw.flatten()
    y_raw = y_raw.flatten()
    z_raw = z_raw.flatten()
    # k_raw = x_raw * y_raw * z_raw
    k_raw = x_raw + y_raw + z_raw

    plot_3d_dataset_as_cube(x_raw, y_raw, z_raw, k_raw,
                            substance_titles=('Alcohol', 'HBr', 'Temperature'),
                            colorbar_title='Yield',
                            npoints=30, sparse_npoints=4, rbf_epsilon=0.04,
                            rbf_smooth=0.001)