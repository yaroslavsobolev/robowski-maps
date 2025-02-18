import plotly.graph_objects as go
import numpy as np

# X, Y, Z = np.mgrid[-1:1:30j, -1:1:30j, -1:1:30j]
# values = np.sin(np.pi * X) * np.cos(np.pi * Z) * np.sin(np.pi * Y)
nsteps = 40
x = np.linspace(0, 10, nsteps)
y = np.linspace(0, 10, nsteps)
z = np.linspace(0, 10, nsteps)
X, Y, Z = np.meshgrid(x, y, z)

# make a 3D gaussian with a peak at point p0, and sigmas sigma0
def gaussian_3d(x,y,z, p0, sigma0):
    return np.exp(-((x-p0[0])**2/(2*sigma0[0]**2) + (y-p0[1])**2/(2*sigma0[1]**2) + (z-p0[2])**2/(2*sigma0[2]**2)))


values = 95*gaussian_3d(X, Y, Z, [1, 3, 4], [2, 4, 2]) + \
          85*gaussian_3d(X, Y, Z, [8, 7, 8], [4, 4, 4]) + \
          85*gaussian_3d(X, Y, Z, [8, 2, 2], [4, 1, 2])

fig = go.Figure(data=go.Volume(
    x=X.flatten(),
    y=Y.flatten(),
    z=Z.flatten(),
    value=values.flatten(),
    isomin=40,
    isomax=90,
    opacity=1,  # max opacity
    opacityscale='max',
    surface_count=6,
    colorscale='viridis'
))
# fig.update_layout(scene = dict(
#                 xaxis_title='Configuration space of all ligands',
#                 yaxis_title='Configuration space of all substrates',
#                 zaxis_title='Conditions space'))

# add arrows
# x_here = 4.9
# add_arrow(fig, [x_here, 1.39, 1], [x_here, 1.39, 0], color='orange', colorscale='Oranges')
#
# add_arrow(fig, [x_here, 1.39, 0], [x_here + 2.5, 1.39, 0], color='magenta', colorscale='Magenta',
#           sizeref=1, linewidth=20, do_line=False)
#
# add_arrow(fig, [x_here + 2.2, 1.39, 0], [x_here + 0.0, 0.1, 0], color='tomato', colorscale='Reds')

fig.update_layout(scene=dict(xaxis_title='Condition parameter A',
                             yaxis_title='Condition parameter B',
                             zaxis_title='Condition parameter C',
                             )
                  )

fig.update_layout(yaxis={'visible': True, 'showticklabels': False})

# set limits of xyz axes to max and min of X, Y, Z
fig.update_layout(scene=dict(xaxis=dict(range=[np.min(X), np.max(X)]),
                             yaxis=dict(range=[np.min(Y), np.max(Y)]),
                             zaxis=dict(range=[np.min(Z), np.max(Z)]),
                             aspectratio=dict(x=1, y=1, z=1)
                             ))

# fig.write_html("figures/generality_illustration.html")
fig.show()
