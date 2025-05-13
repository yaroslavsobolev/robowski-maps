import matplotlib.pyplot as plt

from visualize_results import *
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import interp1d, interp2d

temps_by_run_name: dict = {'simple-reactions/2023-04-11-run01/': 26,
                            'simple-reactions/2023-04-12-run01/': 36,
                            'simple-reactions/2023-04-14-run01/': 16,
                            'simple-reactions/2023-04-14-run02/': 21,
                            'simple-reactions/2023-04-15-run01/': 31,
                            'simple-reactions/2023-04-15-run02/': 11,
                            'simple-reactions/2023-04-28-run01/': 26}
run_name = 'simple-reactions/2023-04-28-run01/'
temperature = temps_by_run_name[run_name]
df_result = pd.read_csv(data_folder + run_name + f'results/product_concentration.csv')
substrates = ['SN1OH01', 'HBr']
show_mystery = True

figsize = (4.5, 8.2)

dilution_factor = 33.3333
df_calibration = df_result.loc[df_result[substrates[1]] == df_result[substrates[1]].min()]
calib_xs = df_calibration['prod1'].to_numpy()
calib_ys = df_calibration[substrates[0]].to_numpy()
plt.plot(calib_xs, calib_ys, 'o-')
calib_xs = np.insert(calib_xs, 0, 0)
calib_ys = np.insert(calib_ys, 0, 0)
plt.plot(calib_xs[:3], calib_ys[:3], 'o-')
plt.show()
points_to_take = 6
abs_to_conc = interp1d(x=calib_xs[:points_to_take],
                       y=calib_ys[:points_to_take], kind='linear', fill_value='extrapolate')

# df_result['prod1'] = df_result['prod1'].apply(lambda x: 0 if x < 0.02 else x)

xs = df_result[substrates[0]].to_numpy()
ys = df_result[substrates[1]].to_numpy()
zs = df_result['prod1'].to_numpy()

# convert to yields
zs = (xs - abs_to_conc(zs)) / ( xs )
zs[zs<0] = 0
print(f'Maximum yield: {max(zs)}')

interp = LinearNDInterpolator((xs, ys), zs)
X = np.linspace(min(xs), max(xs), 300)
Y = np.linspace(min(ys), max(ys), 600)
X, Y = np.meshgrid(X, Y)
Z = interp(X, Y)

# # nonlinear interpolator that is deprecated in scipy
# X = np.linspace(min(xs), max(xs), 300)
# Y = np.linspace(min(ys), max(ys), 600)
# interp_here = interp2d(xs, ys, zs, kind='quintic')
# Z = interp_here(X, Y)

f1 = plt.figure(figsize=figsize)
plt.pcolormesh(X, Y, Z * 100, shading='auto', vmin=0, vmax=100)
plt.plot(xs, ys, "ok", label="input point")

# CS = plt.contour(X, Y, Z, cmap='Greys_r')
# plt.clabel(CS, inline=True, fontsize=10)

plt.title(f'SN1 yields, measured with 33x dilution, \n Temperature {temperature} C (on oven controls),\n run: {run_name}')
cbar = plt.colorbar()
cbar.ax.set_ylabel('Yield, %', rotation=270)
plt.axis("equal")
plt.xlabel(f'Concentration of {substrates[0]}, mol/L')
plt.ylabel(f'Concentration of {substrates[1]}, mol/L')
plt.tight_layout()
plt.gca().set(xlim=((np.min(xs), np.max(xs))), ylim=((np.min(ys), np.max(ys))))
plt.gcf().savefig(data_folder + run_name + f'results/interpolated_yield.png', dpi=300)
plt.show()



# intermediate concentration plot

xs = df_result[substrates[0]].to_numpy()
ys = df_result[substrates[1]].to_numpy()
zs = df_result['carbocat'].to_numpy()
zs = zs - np.median(zs)

# convert to yields
# zs = (xs * dilution_factor - abs_to_conc(zs)) / ( xs * dilution_factor )
# print(f'Maximum yield: {max(zs)}')

interp = LinearNDInterpolator((xs, ys), zs)
X = np.linspace(min(xs), max(xs), 300)
Y = np.linspace(min(ys), max(ys), 600)
X, Y = np.meshgrid(X, Y)
Z = interp(X, Y)

f2 = plt.figure(figsize=figsize)
plt.pcolormesh(X, Y, Z, shading='auto')
plt.plot(xs, ys, "ok", label="input point")

plt.title(f'Reaction intermediate, measured before dilution,\n Temperature {temperature} C (on oven controls),\n run: {run_name}')
cbar = plt.colorbar()
cbar.ax.set_ylabel('Absorbance of intermediate, absorbance units', rotation=270)
plt.axis("equal")
plt.xlabel(f'Concentration of {substrates[0]}, mol/L')
plt.ylabel(f'Concentration of {substrates[1]}, mol/L')
plt.tight_layout()
plt.gca().set(xlim=((np.min(xs), np.max(xs))), ylim=((np.min(ys), np.max(ys))))
plt.gcf().savefig(data_folder + run_name + f'results/reaction_intermediate.png', dpi=300)
plt.show()

if show_mystery:
    # intermediate concentration plot

    xs = df_result[substrates[0]].to_numpy()
    ys = df_result[substrates[1]].to_numpy()
    zs = df_result['unknown-product'].to_numpy()
    zs = zs - np.median(zs)

    # convert to yields
    # zs = (xs * dilution_factor - abs_to_conc(zs)) / ( xs * dilution_factor )
    # print(f'Maximum yield: {max(zs)}')

    interp = LinearNDInterpolator((xs, ys), zs)
    X = np.linspace(min(xs), max(xs), 300)
    Y = np.linspace(min(ys), max(ys), 600)
    X, Y = np.meshgrid(X, Y)
    Z = interp(X, Y)

    f2 = plt.figure(figsize=figsize)
    plt.pcolormesh(X, Y, Z, shading='auto')
    plt.plot(xs, ys, "ok", label="input point")

    plt.title(
        f'Mysterious product, measured before dilution,\n Temperature {temperature} C (on oven controls),\n run: {run_name}')
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Absorbance of unknown product, absorbance units', rotation=270)
    plt.axis("equal")
    plt.xlabel(f'Concentration of {substrates[0]}, mol/L')
    plt.ylabel(f'Concentration of {substrates[1]}, mol/L')
    plt.tight_layout()
    plt.gca().set(xlim=((np.min(xs), np.max(xs))), ylim=((np.min(ys), np.max(ys))))
    plt.gcf().savefig(data_folder + run_name + f'results/reaction_unknown-product.png', dpi=300)
    plt.show()