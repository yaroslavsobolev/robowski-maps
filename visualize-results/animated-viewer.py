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

contourvalues = [0.05, 0.12]
indices_of_outliers = [909, 917]

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
experiment_name = 'multicomp-reactions/2023-03-29-run01/'
# df_results = pd.read_csv(data_folder + experiment_name + f'results/timepoint{timepoint_id:03d}-reaction_results.csv')

df_results = join_data_from_runs(['multicomp-reactions/2023-03-20-run01/',
                                  'multicomp-reactions/2023-03-29-run01/',
                                  'multicomp-reactions/2023-03-31-run01/',
                                  'multicomp-reactions/2023-04-11-run01/'])
df_results.drop(indices_of_outliers, inplace=True)

substances = ['ic001','am001','ald001','ptsa']
product = 'IIO029A'

substrate_cs = []
for substance in substances:
    substrate_cs.append(df_results[substance].to_numpy())

xs0, ys0, zs0, cats = substrate_cs

print('Max concentrations of substrates: ')
for x in [xs0, ys0, zs0]:
    print(max(x))

minimal_concentration_of_substrates = np.min(np.array([xs0, ys0, zs0]))

unique_cats = sorted(list(set(list(cats))))
print(f'Unique cats: {unique_cats}')

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
for cat_here in unique_cats:
    mask = (cats == cat_here)
    xs = xs0[mask]
    ys = ys0[mask]
    zs = zs0[mask]
    ks = ks0[mask]
    # rbf4 = Rbf(xs, ys, zs, ks, epsilon=0.04, smooth=0.00001)  # function="thin_plate"
    rbf4 = Rbf(xs, ys, zs, ks, epsilon=0.04, smooth=0.00008)  # function="thin_plate"
    ti = np.linspace(0, max_ks0, 10)
    npoints = 30j # variable that controls the resolution of the plot
    xnew, ynew, znew = np.mgrid[minimal_concentration_of_substrates:max_xs0:npoints,
                                minimal_concentration_of_substrates:max_ys0:npoints,
                                minimal_concentration_of_substrates:max_zs0:npoints]
    wnew = rbf4(xnew, ynew, znew)
    wnew[wnew<0] = 0
    all_wnews.append(wnew)
    all_points.append((xs, ys, zs, ks))


def curve(n_mer):
    """
    Return the x, y, z coordinates of the points for frame n_mer.
    Parameters
    ----------
    n_mer: int
        Frame number. This is the index of the catalyst concentration in the unique_cats list.

    Returns
    -------

    """
    return all_wnews[n_mer], all_points[n_mer]


class MyModel(HasTraits):
    pTSA_concentration_id = Range(0, len(unique_cats)-1, 17, ) # slider for catalyst concentration
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
            self.plot = self.scene.mlab.contour3d(xnew, ynew, znew, wnew, extent=[0.12, 0.30, 0.12, 0.30, 0.12, 0.30],
                                                  contours=contourvalues, opacity=0.5, vmin=0, vmax=max_ks0,
                                  colormap='summer')
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
            self.scene.mlab.xlabel(f'{substances[0]}')
            self.scene.mlab.ylabel(f'{substances[1]}')
            self.scene.mlab.zlabel(f'{substances[2]}')
            self.scene.mlab.outline(self.plot)
            self.texthere = self.scene.mlab.text3d(max_xs0 * 0.7, max_ys0*1.3, max_zs0/2, 'Catalyst', scale=0.005)
            self.vslice = self.scene.mlab.volume_slice(xnew, ynew, znew, wnew, plane_orientation='x_axes', opacity=0.5,
                                                       vmin=0, vmax=max_ks0)#, colormap='summer')
            self.scene.background = (1, 1, 1)
            cb = self.scene.mlab.colorbar(object=self.vslice, title="Yield")
            cb.scalar_bar.unconstrained_font_size = True
            cb.label_text_property.font_size = 15
            ax1.axes.font_factor = 0.83
            xs, ys, zs, ks = the_points
            self.plot_points = self.scene.mlab.points3d(xs, ys, zs, ks, vmin=0, vmax=max_ks0, resolution=16)
        else:
            # Explanation for the following code:
            self.plot.mlab_source.trait_set(x=xnew, y=ynew, z=znew, scalars=wnew)
            self.plot.contour.contours = contourvalues
            self.vslice.mlab_source.trait_set(x=xnew, y=ynew, z=znew, scalars=wnew)
            xs, ys, zs, ks = the_points
            self.plot_points.mlab_source.reset(x=xs, y=ys, z=zs, scalars=ks)
            self.texthere.text = f'Catalyst (pTSA) {unique_cats[self.pTSA_concentration_id]:.3f} mol/L'
            self.texthere.vector_text.update()
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