import plotly.graph_objects as go
import numpy as np
import matplotlib.pyplot as plt

# nsteps = 400
# x = np.arange(0, 30, 1)
#
# def gaussian_1d_curve(x, p0, sigma0):
#     return np.exp(-((x - p0) ** 2 / (2 * sigma0 ** 2)))
#
# y = 100*gaussian_1d_curve(x, 5, 2) + 30*gaussian_1d_curve(x, 22, 4) + 30*gaussian_1d_curve(x, 30, 6)
#
# # make bar plot vertical
# plt.bar(x, y, width=1, color='grey', edgecolor='black')
#
# plt.xlabel('Conditions')
# plt.ylabel('Yield, %')
# # make ticks every 1 point
# plt.xlim(-0.5, 29.5)
# plt.show()

nsteps = 400
x = np.arange(0, 15, 1)
y = np.arange(0, 10, 1)

def gaussian_2d_surface(x, y, p0, sigma0):
    return np.exp(-((x - p0[0]) ** 2 / (2 * sigma0[0] ** 2) + (y - p0[1]) ** 2 / (2 * sigma0[1] ** 2)))

# plot and make plot with pcolor, nearest neighbour
X, Y = np.meshgrid(x, y)
Z = 100*gaussian_2d_surface(X, Y, [0, 1], [2, 2]) + 90*gaussian_2d_surface(X, Y, [8, 6], [4, 2]) + 80*gaussian_2d_surface(X, Y, [14, 0], [3, 2])
fig, ax = plt.subplots()
c = ax.pcolormesh(X, Y, Z, cmap='viridis', vmin=0, vmax=100, shading='nearest')
plt.xlabel('Condition parameter A')
plt.ylabel('Condition parameter B')
cbar = fig.colorbar(c, label='Yield, %')
plt.show()