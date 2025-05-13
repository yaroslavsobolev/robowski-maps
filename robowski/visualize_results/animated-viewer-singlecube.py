from traits.api import HasTraits, Range, Instance, \
        on_trait_change
from traitsui.api import View, Item, Group
from mayavi.core.api import PipelineBase
from mayavi.core.ui.api import MayaviScene, SceneEditor, \
                MlabSceneModel
from tvtk.tools import visual
from scipy.interpolate import Rbf
from visualize_results import *
import time
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import interp1d

os.environ['ETS_TOOLKIT'] = 'qt'
os.environ['QT_API'] = 'pyqt5'

organize_run_results = importlib.import_module("misc_scripts.organize_run_results")

contourvalues = [0.5, 0.7, 0.9]
indices_of_outliers = []

def Arrow_From_A_to_B(x1, y1, z1, x2, y2, z2):
    """
    Create an arrow from A to B given by coordinates x1, y1, z1 for A and x2, y2, z2 for B
    """
    ar1=visual.arrow(x=x1, y=y1, z=z1)
    ar1.color = (0.5, 0.5, 0.5)

    arrow_length=np.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
    ar1.length_cone=0.12*30/arrow_length
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

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

# timepoint_id = 1
experiment_name = 'simple-reactions/2023-04-15-run02/'
# df_results = pd.read_csv(data_folder + experiment_name + f'results/timepoint{timepoint_id:03d}-reaction_results.csv')

df_results = organize_run_results.join_data_from_runs(['simple-reactions/2023-04-11-run01/',
                            'simple-reactions/2023-04-12-run01/',
                            'simple-reactions/2023-04-14-run01/',
                            'simple-reactions/2023-04-14-run02/',
                            'simple-reactions/2023-04-15-run01/',
                            'simple-reactions/2023-04-15-run02/'],
                                 round_on_columns=('SN1OH01', 'HBr', 'Temperature'))
df_results.drop(indices_of_outliers, inplace=True)

substances = ['SN1OH01', 'HBr', 'Temperature']
substance_titles = ['Alcohol', 'HBr', 'Temperature']
substrates = ['SN1OH01', 'HBr']
product = 'prod1'

yscale = 0.33
dilution_factor = 33.3333
df_calibration = df_results.loc[(df_results[substrates[1]] == df_results[substrates[1]].min()) &
                                (df_results['Temperature'] == df_results['Temperature'].min())]
calib_xs = df_calibration['prod1'].to_numpy()
calib_ys = df_calibration[substrates[0]].to_numpy()
plt.plot(calib_xs, calib_ys, 'o-')
calib_xs = np.insert(calib_xs, 0, 0)
calib_ys = np.insert(calib_ys, 0, 0)
plt.plot(calib_xs[:3], calib_ys[:3], 'o-')
plt.show()
abs_to_conc = interp1d(x=calib_xs[:],
                       y=calib_ys[:], kind='linear', fill_value='extrapolate')

# df_result['prod1'] = df_result['prod1'].apply(lambda x: 0 if x < 0.02 else x)

xs = df_results[substrates[0]].to_numpy()
ys = df_results[substrates[1]].to_numpy()
zs = df_results['prod1'].to_numpy()

# convert to yields
zs = (xs - abs_to_conc(zs)) / ( xs )
zs[zs<0] = 0
df_results['yield'] = zs
print(f'Maximum yield: {max(zs)}')

# convert to mM

for substrate in substrates:
    df_results[substrate] = df_results[substrate].apply(lambda x: x * 1000)
    df_results[substrate] = df_results[substrate].apply(lambda x: x if x>1e-10 else 0)
df_results['carbocat'] = df_results['carbocat'].apply(lambda x: x * 8)
df_results[substrates[1]] = df_results[substrates[1]].apply(lambda x: x * yscale)

substrate_cs = []
for substance in substances:
    substrate_cs.append(df_results[substance].to_numpy())

xs0, ys0, zs0 = substrate_cs

print('Max concentrations of substrates: ')
for x in [xs0, ys0, zs0]:
    print(max(x))

minimal_concentration_of_substrates = np.min(np.array([xs0, ys0, zs0]))


yields = df_results['yield'].to_numpy()

print('Min and max yields:')
ks0 = yields
print(max(ks0))
print(min(ks0))

max_ks0 = np.max(ks0)
max_xs0 = np.max(xs0)
max_ys0 = np.max(ys0)
max_zs0 = np.max(zs0)

# We have 3 substrate concentrations and 1 catalyst concentration. We want to plot the yield in
# 3D, so we need to interpolate the data. We do that by creating a RBF function for each catalyst
# concentration, and then plotting the results. So we first create a list of all the RBF functions,
# and then plot them based on the current value of the catalyst concentration slider.
# The results of RBF interpolation will be stored in all_wnews variable.
all_wnews = []
all_points = []
unique_cats = ['yield', 'carbocat', 'unknown-product']

npoints = 30j  # variable that controls the resolution of the plot
xnew, ynew, znew = np.mgrid[minimal_concentration_of_substrates:max_xs0:npoints,
                   minimal_concentration_of_substrates:max_ys0:npoints,
                   minimal_concentration_of_substrates:max_zs0:npoints]

# construct sparse points for spheres
sparse_npoints = [4j, 4j, 4j]
xnew_sparse, ynew_sparse, znew_sparse = np.mgrid[np.min(xs0):max_xs0:sparse_npoints[0],
                   np.min(ys0):max_ys0:sparse_npoints[1],
                   np.min(zs0):max_zs0:sparse_npoints[2]]
# sparse_mask = (xnew_sparse == max_xs0) | \
#               (ynew_sparse == minimal_concentration_of_substrates) | \
#               (znew_sparse == minimal_concentration_of_substrates)
# xnew_sparse = xnew_sparse[sparse_mask]
# ynew_sparse = ynew_sparse[sparse_mask]
# znew_sparse = znew_sparse[sparse_mask]

for cat_here in unique_cats:
    xs = xs0
    ys = ys0
    zs = zs0
    ks = df_results[cat_here].to_numpy()

    # RBF version
    # rbf4 = Rbf(xs, ys, zs, ks, epsilon=0.04, smooth=0.00001)  # function="thin_plate"
    # rbf4 = Rbf(xs, ys, zs, ks, epsilon=0.04, smooth=0.00008)  # function="thin_plate"
    rbf4 = Rbf(xs, ys, zs, ks, epsilon=3, smooth=5)  # function="thin_plate"
    wnew = rbf4(xnew, ynew, znew)

    # # Linead ND version
    interp_here = LinearNDInterpolator((xs, ys, zs), ks)
    wnew = interp_here((xnew, ynew, znew))

    wnew[wnew<0] = 0
    all_wnews.append(wnew)

    # THis uses the raw points from dataframe
    sparse_ks = interp_here((xnew_sparse, ynew_sparse, znew_sparse))
    # all_points.append((xs, ys, zs, ks))
    all_points.append((xnew_sparse, ynew_sparse, znew_sparse, sparse_ks))

    # This uses the points taken from RBF already


def curve(n_mer):
    return all_wnews[n_mer], all_points[n_mer]


class MyModel(HasTraits):
    pTSA_concentration_id = Range(0, len(unique_cats)-1, 0, ) # slider for catalyst concentration
    second_slider_id = Range(0, 3, )  # slider for catalyst concentration
    scene = Instance(MlabSceneModel, ())
    scene.background = (1, 1, 1)
    plot = Instance(PipelineBase)

    # When the scene is activated, or when the parameters are changed, we
    # update the plot.
    @on_trait_change('second_slider_id,scene.activated')
    def do_animation(self):
        if self.second_slider_id == 2:
            for i in range(len(unique_cats)):
                time.sleep(1)
                self.pTSA_concentration_id = i
                self.update_plot()


    @on_trait_change('pTSA_concentration_id,scene.activated')
    def update_plot(self):
        wnew, the_points = curve(self.pTSA_concentration_id)
        if self.plot is None:
            self.plot = self.scene.mlab.contour3d(xnew, ynew, znew, wnew, extent=[np.min(xs0), max_xs0,
                                                                                  np.min(ys0), max_ys0,
                                                                                  np.min(zs0), max_zs0],
                                                  contours=contourvalues, opacity=1, vmin=0, vmax=max_ks0)  # ,
            # colormap='summer')
            self.plot.actor.actor.property.ambient = 0.0
            self.fig = self.scene.mlab.gcf()
            visual.set_viewer(self.fig)
            for i in range(3):
                start = np.array([np.min(xs0), np.min(ys0), np.min(zs0)])
                end = np.array([np.min(xs0), np.min(ys0), np.min(zs0)])
                end[i] = list([max_xs0, max_ys0, max_zs0])[i]
                self.arr = Arrow_From_A_to_B(start[0], start[1], start[2], end[0], end[1], end[2])
            self.arr_temp = Arrow_From_A_to_B(np.max(xs0), np.min(ys0), np.min(zs0),
                                              np.max(xs0), np.min(ys0), np.max(zs0))
            ax1 = self.scene.mlab.axes(color=(0.5, 0.5, 0.5), nb_labels=0, ranges=[np.min(xs0), max_xs0,
                                                                                  np.min(ys0)/yscale, max_ys0/yscale,
                                                                                  np.min(zs0), max_zs0],
                                       x_axis_visibility=False, y_axis_visibility=False, z_axis_visibility=False)
            self.scene.foreground = (0.2, 0.2, 0.2)
            self.scene.mlab.xlabel(f'{substance_titles[0]}')
            self.scene.mlab.ylabel(f'{substance_titles[1]}')
            self.scene.mlab.zlabel(f'{substance_titles[2]}')
            self.scene.mlab.outline(self.plot)
            self.texthere = self.scene.mlab.text3d(max_xs0 * (0.2), max_ys0 * 0.6, max_zs0 / 2, 'Catalyst', scale=0.005)
            # self.vslice = self.scene.mlab.volume_slice(xnew, ynew, znew, wnew, plane_orientation='x_axes',
            #                                            plane_opacity=0.01, opacity=0.01, transparent=True, #opacity=0.5,
            #                                            vmin=0, vmax=max_ks0)#, colormap='summer')
            self.scene.background = (1, 1, 1)
            cb = self.scene.mlab.colorbar(object=self.plot, title="Yield")
            cb.scalar_bar.unconstrained_font_size = True
            cb.label_text_property.font_size = 15
            ax1.axes.font_factor = 0.83
            xs, ys, zs, ks = the_points
            self.plot_points = self.scene.mlab.points3d(xs, ys, zs, ks, vmin=0, vmax=max_ks0, resolution=16, scale_factor=4)  # ,
            # colormap='summer')

            # self.scene.scene.camera.position = [-0.168916619662477, 0.3818710137155461, 0.8890913200020035]
            # self.scene.scene.camera.focal_point = [0.21530580531460342, 0.21442896061186448, 0.21702529116995473]
            # self.scene.scene.camera.view_angle = 30.0
            # self.scene.scene.camera.view_up = [-0.8744518469052971, -0.11364676850842695, -0.47161253105860845]
            # self.scene.scene.camera.clipping_range = [0.42598520915716326, 1.257159582158653]
            # self.scene.scene.camera.compute_view_plane_normal()
            # self.scene.scene.render()

        else:
            # Explanation for the following code:
            # self.plot.mlab_source.trait_set(x=xnew, y=ynew, z=znew, scalars=wnew)
            # self.plot.contour.contours = contourvalues
            self.plot.remove()
            self.plot = self.scene.mlab.contour3d(xnew, ynew, znew, wnew, extent=[np.min(xs0), max_xs0,
                                                                                  np.min(ys0), max_ys0,
                                                                                  np.min(zs0), max_zs0],
                                                  contours=contourvalues, opacity=1, vmin=0, vmax=max_ks0)  # ,
            # colormap='summer')
            self.plot.actor.actor.property.ambient = 0.0

            # self.vslice.mlab_source.trait_set(x=xnew, y=ynew, z=znew, scalars=wnew, plane_opacity=0.01, opacity=0.01, transparent=True)
            xs, ys, zs, ks = the_points
            self.plot_points.mlab_source.reset(x=xs, y=ys, z=zs, scalars=ks)
            self.texthere.text = f'Catalyst (pTSA) {unique_cats[self.pTSA_concentration_id]:.3f} mol/L'
            self.texthere.vector_text.update()
            # # update the camera positino
            # pos, fp, vu = cam_pos_and_focalpoint_by_blendfactor(blend_factor=(self.pTSA_concentration_id / 500) ** 2)
            # self.scene.scene.camera.position = pos
            # self.scene.scene.camera.focal_point = fp
            # self.scene.scene.camera.view_up = vu
            # self.scene.scene.render()

            create_folder_unless_it_exists(data_folder + experiment_name + f'results/4d-viewer-frames')
            self.scene.mlab.savefig(
                data_folder + experiment_name + f'results/4d-viewer-frames/{self.pTSA_concentration_id:05d}.png')
            # time.sleep(1)
            # self.pTSA_concentration_id += 1
            # self.update_plot()


    # The layout of the dialog created
    # Explanation for the following code:
    view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                     height=250, width=300, show_label=False),
                Group(
                        '_', 'pTSA_concentration_id', 'second_slider_id',
                     ),
                resizable=True,
                )

# Explanation for the following code:
my_model = MyModel()
my_model.configure_traits()