import pickle

from traits.api import HasTraits, Range, Instance, \
        on_trait_change
from traitsui.api import View, Item, Group
from mayavi.core.api import PipelineBase
from mayavi.core.ui.api import MayaviScene, SceneEditor, \
                MlabSceneModel
from tvtk.tools import visual
from scipy.interpolate import Rbf, LinearNDInterpolator
from visualize_results import *

os.environ['ETS_TOOLKIT'] = 'qt'
os.environ['QT_API'] = 'pyqt5'

organize_run_results = importlib.import_module("misc_scripts.organize_run_results")

contourvalues = [0.10, 0.14]


pos1 = [-0.168916619662477, 0.3818710137155461, 0.8890913200020035]
fp1 = [0.21530580531460342, 0.21442896061186448, 0.21702529116995473]
vu1 = [-0.8744518469052971, -0.11364676850842695, -0.47161253105860845]

pos2 = [-0.1297603269375693, 0.9093667836146796, 0.37617065442515396]
fp2 = [0.21530580531460342, 0.21442896061186448, 0.21702529116995473]
vu2 = [-0.9000865671921399, -0.4230780617551936, -0.10415913412532742]

def cam_pos_and_focalpoint_by_blendfactor(blend_factor):
    pos = np.array(pos1) * (1 - blend_factor) + np.array(pos2) * blend_factor
    fp = np.array(fp1) * (1 - blend_factor) + np.array(fp2) * blend_factor
    vu = np.array(vu1) * (1 - blend_factor) + np.array(vu2) * blend_factor
    return pos.tolist(), fp.tolist(), vu.tolist()

def Arrow_From_A_to_B(x1, y1, z1, x2, y2, z2):
    """
    Create an arrow from A to B given by coordinates x1, y1, z1 for A and x2, y2, z2 for B
    """
    ar1=visual.arrow(x=x1, y=y1, z=z1)
    ar1.length_cone=0.12
    ar1.color = (0.5, 0.5, 0.5)

    arrow_length=np.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
    ar1.radius_shaft = 0.015
    ar1.radius_cone = 0.04
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
experiment_name = 'multicomp-reactions/2023-06-19-run01/'
# df_results = pd.read_csv(data_folder + experiment_name + f'results/timepoint{timepoint_id:03d}-reaction_results.csv')

# list_of_runs = tuple(['2023-06-20-run01',
#                     '2023-06-21-run01',
#                     '2023-06-21-run02',
#                     '2023-06-22-run01',
#                     '2023-06-22-run02',
#                     '2023-06-22-run03',
#                     '2023-06-23-run01',
#                     '2023-06-23-run02',
#                     '2023-06-26-run01',
#                     '2023-06-26-run02',
#                     '2023-06-27-run01',
#                     '2023-06-27-run02',
#                     '2023-06-27-run03',
#                     '2023-06-28-run01',
#                     '2023-06-28-run02',
#                     '2023-06-28-run03'])
#
# df_results = organize_run_results.join_data_from_runs([f'multicomp-reactions/{run}/' for run in list_of_runs])

df_results = pd.read_csv(data_folder + experiment_name + f'results/product_concentration_after_substituting_outliers.csv')

print(f"There are {df_results[df_results['is_outlier'] == 1].shape[0]} outliers.")
df_results = df_results[df_results['is_outlier'] == 0]


# df_results = pd.read_csv(data_folder + experiment_name + f'results/interpolated_product_concentration.csv')
# replace negative values in yield column with zeros
df_results['yield'] = df_results['yield'].apply(lambda x: 0 if x < 0 else x)
# df_results.drop('Unnamed: 0', inplace=True, axis=1)

substances = ['ic001','am001','ald001','ptsa']
substance_titles = ['Isocyanide', 'Amine', 'Aldehyde', 'p-TSA']
# substance_titles = ['', '', '', '']

padding_rows_count = (df_results[substances] == 0).all(axis=1).sum()
print(f"There are {padding_rows_count} padding rows (with zero concentrations of substrates).")
df_results = df_results[(df_results[substances] != 0).any(axis=1)]

substrate_cs = []
for substance in substances:
    substrate_cs.append(df_results[substance].to_numpy())

xs0, ys0, zs0, cats = substrate_cs

print('Max concentrations of substrates: ')
for x in [xs0, ys0, zs0]:
    print(max(x))

minimal_concentration_of_substrates = np.min(np.array([xs0, ys0, zs0]))

unique_cats = list(sorted(list(df_results['ptsa'].unique())))
# pickle unique_cats to file
with open(data_folder + experiment_name + 'results/unique_cats.pickle', 'wb') as f:
    pickle.dump(unique_cats, f)

print(f'{len(unique_cats)} unique cats: {unique_cats}')

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

npoints = 30j  # variable that controls the resolution of the plot
xnew, ynew, znew = np.mgrid[minimal_concentration_of_substrates:max_xs0:npoints,
                   minimal_concentration_of_substrates:max_ys0:npoints,
                   minimal_concentration_of_substrates:max_zs0:npoints]

# construct sparse points for spheres
sparse_npoints = 8j
xnew_sparse, ynew_sparse, znew_sparse = np.mgrid[minimal_concentration_of_substrates:max_xs0:sparse_npoints,
                   minimal_concentration_of_substrates:max_ys0:sparse_npoints,
                   minimal_concentration_of_substrates:max_zs0:sparse_npoints]
sparse_mask = (xnew_sparse == max_xs0) | \
              (ynew_sparse == minimal_concentration_of_substrates) | \
              (znew_sparse == minimal_concentration_of_substrates)
xnew_sparse = xnew_sparse[sparse_mask]
ynew_sparse = ynew_sparse[sparse_mask]
znew_sparse = znew_sparse[sparse_mask]

for cat_here in unique_cats:
    mask = (cats == cat_here)
    xs = xs0[mask]
    ys = ys0[mask]
    zs = zs0[mask]
    ks = ks0[mask]

    # RBF version
    # rbf4 = Rbf(xs, ys, zs, ks, epsilon=0.04, smooth=0.00001)  # function="thin_plate"
    # rbf4 = Rbf(xs, ys, zs, ks, epsilon=0.04, smooth=0.00008)  # function="thin_plate"
    rbf4 = Rbf(xs, ys, zs, ks, epsilon=0.04, smooth=0.8)  # function="thin_plate"
    wnew = rbf4(xnew, ynew, znew)

    # # Linead ND version
    # interp_here = LinearNDInterpolator((xs, ys, zs), ks)
    # wnew = interp_here((xnew, ynew, znew))

    wnew[wnew<0] = 0
    all_wnews.append(wnew)

    # # THis uses the raw points from dataframe
    all_points.append((xs, ys, zs, ks))

    # this uses sparsely resampled points
    # sparse_ks = rbf4(xnew_sparse, ynew_sparse, znew_sparse)
    # all_points.append((xnew_sparse, ynew_sparse, znew_sparse, sparse_ks))

def curve(n_mer):
    return all_wnews[n_mer], all_points[n_mer]


class MyModel(HasTraits):
    pTSA_concentration_id = Range(0, len(unique_cats)-1, 1) # slider for catalyst concentration
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
                # time.sleep(1)
                self.pTSA_concentration_id = i
                self.update_plot()


    @on_trait_change('pTSA_concentration_id,scene.activated')
    def update_plot(self):
        wnew, the_points = curve(self.pTSA_concentration_id)
        if self.plot is None:
            self.plot = self.scene.mlab.contour3d(xnew, ynew, znew, wnew, extent=[0.12, 0.30, 0.12, 0.30, 0.12, 0.30],
                                                  contours=contourvalues, opacity=1, vmin=0, vmax=max_ks0)#,
                                  # colormap='summer')
            self.plot.actor.actor.property.ambient = 0.5
            self.fig = self.scene.mlab.gcf()
            visual.set_viewer(self.fig)
            for i in range(3):
                start = [minimal_concentration_of_substrates]*3
                end = [minimal_concentration_of_substrates]*3
                end[i] = list([max_xs0, max_ys0, max_zs0])[i]
                self.arr = Arrow_From_A_to_B(start[0], start[1], start[2], end[0], end[1], end[2])
            ax1 = self.scene.mlab.axes(color=(0.5, 0.5, 0.5), nb_labels=4)
            self.scene.foreground = (0.2, 0.2, 0.2)
            self.scene.mlab.xlabel(f'{substance_titles[0]}')
            self.scene.mlab.ylabel(f'{substance_titles[1]}')
            self.scene.mlab.zlabel(f'{substance_titles[2]}')
            self.scene.mlab.outline(self.plot)
            self.texthere = self.scene.mlab.text3d(max_xs0 * (0.2), max_ys0 * 0.6, max_zs0/2*1.3, 'Catalyst', scale=0.005)
            # self.vslice = self.scene.mlab.volume_slice(xnew, ynew, znew, wnew, plane_orientation='x_axes',
            #                                            plane_opacity=0.01, opacity=0.01, transparent=True, #opacity=0.5,
            #                                            vmin=0, vmax=max_ks0)#, colormap='summer')
            self.scene.background = (1, 1, 1)
            cb = self.scene.mlab.colorbar(object=self.plot, title="Yield", label_fmt='%.2f', nb_labels=4)
            cb.scalar_bar.unconstrained_font_size = True
            cb.label_text_property.font_size = 15
            ax1.axes.font_factor = 0.83
            xs, ys, zs, ks = the_points
            self.plot_points = self.scene.mlab.points3d(xs, ys, zs, ks, vmin=0, vmax=max_ks0, resolution=16)#,
                                                        # colormap='summer')

            self.scene.scene.camera.position = [-0.168916619662477, 0.3818710137155461, 0.8890913200020035]
            self.scene.scene.camera.focal_point = [0.21530580531460342, 0.21442896061186448, 0.21702529116995473]
            self.scene.scene.camera.view_angle = 30.0
            self.scene.scene.camera.view_up = [-0.8744518469052971, -0.11364676850842695, -0.47161253105860845]
            self.scene.scene.camera.clipping_range = [0.42598520915716326, 1.257159582158653]
            self.scene.scene.camera.compute_view_plane_normal()
            self.scene.scene.render()

        else:
            # self.plot.mlab_source.trait_set(x=xnew, y=ynew, z=znew, scalars=wnew)
            # self.plot.contour.contours = contourvalues
            # self.plot.remove()
            # self.plot = self.scene.mlab.contour3d(xnew, ynew, znew, wnew, extent=[0.12, 0.30, 0.12, 0.30, 0.12, 0.30],
            #                                       contours=contourvalues, opacity=1, vmin=0, vmax=max_ks0)
                                  # colormap='summer')
            self.plot.actor.actor.property.ambient = 0.2
            self.plot.actor.actor.property.specular = 0.2


            # self.vslice.mlab_source.trait_set(x=xnew, y=ynew, z=znew, scalars=wnew, plane_opacity=0.01, opacity=0.01, transparent=True)
            xs, ys, zs, ks = the_points
            self.plot_points.mlab_source.reset(x=xs, y=ys, z=zs, scalars=ks)
            self.texthere.text = f'Catalyst (p-Toluenesulfonic acid): {unique_cats[self.pTSA_concentration_id]:.3f} mol/L'
            self.texthere.vector_text.update()
            # update the camera positino
            # pos, fp, vu = cam_pos_and_focalpoint_by_blendfactor(blend_factor=(self.pTSA_concentration_id/500)**2)
            pos, fp, vu = cam_pos_and_focalpoint_by_blendfactor(blend_factor=(self.pTSA_concentration_id / 30) ** 2)
            self.scene.scene.camera.position = pos
            self.scene.scene.camera.focal_point = fp
            self.scene.scene.camera.view_up = vu
            self.scene.scene.render()

            create_folder_unless_it_exists(data_folder + experiment_name + f'results/4d-viewer-frames')
            self.scene.mlab.savefig(data_folder + experiment_name + f'results/4d-viewer-frames/{self.pTSA_concentration_id:05d}.png')
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