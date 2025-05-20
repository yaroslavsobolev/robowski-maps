from robowski.settings import *
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib.ticker import FuncFormatter
from visualize_results import *
from scipy.interpolate import Rbf, LinearNDInterpolator

# import 3d matplotlib
from mpl_toolkits.mplot3d import Axes3D

# df_results = join_data_from_runs(['multicomp-reactions/2023-03-20-run01/',
#                                   'multicomp-reactions/2023-03-29-run01/',
#                                   'multicomp-reactions/2023-03-31-run01/',
#                                   'multicomp-reactions/2023-04-11-run01/'])
run_name = 'multicomp-reactions/2023-06-19-run01/'

# # dataset with outliers remeasured
# df_results = pd.read_csv(data_folder + run_name + f'results/product_concentration_after_substituting_outliers.csv')

# data with additional column containing the values predicted by the best-fit model (v3)
df_results = pd.read_csv(data_folder + run_name + f'results/ugi_v3_model_yields_for_entire_exp_dataset_2025-01-31.csv')

# # model results for model after two iterations of optimization on the entire raw dataset
# df_results_fullopt = pd.read_csv(data_folder + run_name + f'results/ugi_v3_model_yields_for_entire_exp_dataset_2025-02-14.csv')
# df_results['v3_full_2iter_model_yield'] = df_results_fullopt['v3_full_2iter_model_yield']

# plot max yield over the entire dataframe
print(f'Max yield: {df_results["yield"].max()}')

# group dataframe by unique values of 'ptsa' column and count the number of rows in each group
ptsa_unique = df_results['ptsa'].unique()
ptsa_interesting = ptsa_unique[0]
df_interesting = df_results[df_results['ptsa'] == ptsa_interesting]
xs = df_interesting['ic001'].to_numpy()
ys = df_interesting['am001'].to_numpy()
zs = df_interesting['ald001'].to_numpy()

# make 3d plot with the xs, ys, zs
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(xs, ys, zs, c='r', marker='o')
ax.set_xlabel('Isocyanide')
ax.set_ylabel('Amine')
ax.set_zlabel('Aldehyde')
plt.show()

do_plot = False
plot_interpolants = True
substances = ['ic001', 'am001', 'ald001', 'ptsa']
product = 'IIO029A'
catalyst_name = 'ptsa'
relative_std = 0.146

# here I drop two first columns from df_results pandas dataframe because they are not needed
# for the interpolation

# indices_of_outliers = [901, 909, 917]
# df_results.drop(indices_of_outliers, inplace=True)

target_folder = data_folder + run_name + f'results/interpolation_plots/'
# if folder does not exist, create it
if not os.path.exists(target_folder):
    os.makedirs(target_folder)

substrate_cs = []
for substance in substances:
    substrate_cs.append(df_results[substance].to_numpy())

xs0, ys0, zs0, cats = substrate_cs
ks0 = df_results['yield'].to_numpy()

for substance in substances:
    print(f'{substance} min: {np.min(df_results[substance].to_numpy())}')
    print(f'{substance} max: {np.max(df_results[substance].to_numpy())}')
    print(f'{substance} unique: {sorted(df_results[substance].unique())}')

minimal_concentration_of_substrates = np.min(np.array([xs0, ys0, zs0]))

unique_cats = sorted(list(set(list(cats))))
print(f'{len(unique_cats)} unique cats: {unique_cats}')

npoints_over_catalyst_range = 100
new_cats = np.linspace(df_results[catalyst_name].min(), df_results[catalyst_name].max(), npoints_over_catalyst_range)

df_combinations = df_results.groupby(['ic001', 'am001', 'ald001']).size().reset_index().rename(columns={0: 'count'})

df_interpolated = df_results.drop(df_results.index)
df_interpolated = df_interpolated[substances + ['yield'] + ['v3_model_yield']]
for index, row in df_combinations.iterrows():
    print(f'Index: {index}')
    print(row)
    # find all the rows from df_results that have this combination of concentrations
    df_this_combination = df_results[(df_results['ic001'] == row['ic001']) &
                                     (df_results['am001'] == row['am001']) &
                                     (df_results['ald001'] == row['ald001'])]
    # if there is only one row in df_this_combination, then there is no need to interpolate
    if len(df_this_combination) == 1:
        df_temporary = df_interpolated.drop(df_interpolated.index)
        dictionary_here = {substance: row[substance] for substance in substances}
        dictionary_here['yield'] = df_this_combination['yield'].to_numpy()[0]
        df_temporary.loc[0] = dictionary_here
        # df_temporary.loc[0] = [row['ald001'], df_this_combination['ptsa'].to_numpy()[0], row['am001'], row['ic001'],
        #                        df_this_combination['yield'].to_numpy()[0]]
        df_interpolated = df_interpolated.append(df_temporary, ignore_index=True)
        continue

    # if index == 18:
    #     df_this_combination['yield'] = df_this_combination['yield'] * 18.6/10.7
    #
    # if index == 29:
    #     df_this_combination['yield'] = df_this_combination['yield'] * 12.5/8.8

    catalyst_concentration_limits_here = df_this_combination['ptsa'].min(), df_this_combination['ptsa'].max()

    # For catalyst concentrations where this combination was not measured,
    # make RBF interpolation from nearby points in substrate space
    if do_plot:
        plt.errorbar(df_this_combination['ptsa'].to_numpy(), df_this_combination['yield'].to_numpy(),
                     yerr=df_this_combination['yield'].to_numpy() * relative_std, linestyle='None',
                     marker='o', capsize=3, label='data', alpha=0.5, color='C0')
    xnew = row['ic001']
    ynew = row['am001']
    znew = row['ald001']
    for catalyst_concentration in unique_cats:
        if catalyst_concentration < catalyst_concentration_limits_here[0] or catalyst_concentration > \
                catalyst_concentration_limits_here[1]:
            mask = (cats == catalyst_concentration)
            xs = xs0[mask]
            ys = ys0[mask]
            zs = zs0[mask]
            ks = ks0[mask]

            # # RBF version
            # rbf4 = Rbf(xs, ys, zs, ks, epsilon=0.04, smooth=0.8)  # function="thin_plate"
            # wnew = rbf4(xnew, ynew, znew)

            # Linead ND version
            interp_here = LinearNDInterpolator((xs, ys, zs), ks)
            wnew = float(interp_here((xnew, ynew, znew)))
            if wnew < 0:
                wnew = 0
            # append one row to df_this_combination with wnew as yield and catalyst_concentration as ptsa
            df_temporary = df_interpolated.drop(df_interpolated.index)
            dictionary_for_this_row = {substance: row[substance] for substance in substances if
                                       not (substance == catalyst_name)}
            dictionary_for_this_row['yield'] = wnew
            dictionary_for_this_row[catalyst_name] = catalyst_concentration
            df_temporary.loc[0] = dictionary_for_this_row
            df_this_combination = df_this_combination.append(df_temporary, ignore_index=True)
    assert df_this_combination[catalyst_name].min() == unique_cats[0]
    assert df_this_combination[catalyst_name].max() == unique_cats[-1]

    catalyst_concentration_limits_here = df_this_combination['ptsa'].min(), df_this_combination['ptsa'].max()

    # make interpolator of yields over catalyst concentrations
    xs = df_this_combination['ptsa'].to_numpy()
    ys = df_this_combination['yield'].to_numpy()
    yerr = df_this_combination['yield'].to_numpy() * relative_std

    f1 = plt.figure(1, figsize=(5, 4))
    plt.errorbar(xs, ys, yerr=yerr, linestyle='None', marker='x', capsize=3, label='data', alpha=0.5, color='black')
    # tck = make_spline(xs, ys, yerr, do_plot=do_plot, spline_smoothing_factor=0.05)
    tck = make_spline(xs, ys, yerr, do_plot=(plot_interpolants), spline_smoothing_factor=0.0005,
                      layout_threshold=0.02)
    if plot_interpolants:
        xs_new = np.linspace(np.min(xs), np.max(xs), 100)
        ys_new = splev(xs_new, tck)
        plt.plot(xs_new, ys_new, color='black', linewidth=5, alpha=0.3)

    # plot the 'v3_model_yield' column from df_results
    xs_model = df_this_combination['ptsa'].to_numpy()
    ys_model = df_this_combination['v3_model_yield'].to_numpy()
    # sort xs_model and ys_model by increasing values of xs_model
    xs_model, ys_model = zip(*sorted(zip(xs_model, ys_model)))
    plt.plot(xs_model, ys_model, color='C0', linewidth=5, alpha=0.3)

    # xs_model = df_this_combination['ptsa'].to_numpy()
    # ys_model = df_this_combination['v3_full_2iter_model_yield'].to_numpy()
    # # sort xs_model and ys_model by increasing values of xs_model
    # xs_model, ys_model = zip(*sorted(zip(xs_model, ys_model)))
    # plt.plot(xs_model, ys_model, color='C0', linewidth=5, alpha=0.3)

    xs_new = new_cats[np.where(
        (new_cats >= catalyst_concentration_limits_here[0]) & (new_cats <= catalyst_concentration_limits_here[1]))]
    ys_new = splev(xs_new, tck)

    # compute ys_model_new using linear interpolation for xs_new using the x-y data of xs_model and ys_model
    ys_model_new = np.interp(xs_new, xs_model, ys_model)

    df_temporary = df_interpolated.drop(df_interpolated.index)
    for i, ptsa_new in enumerate(xs_new):
        dictionary_for_this_row = {substance: row[substance] for substance in substances if
                                   not (substance == catalyst_name)}
        dictionary_for_this_row['yield'] = ys_new[i]
        dictionary_for_this_row['v3_model_yield'] = ys_model_new[i]
        dictionary_for_this_row[catalyst_name] = ptsa_new
        df_temporary.loc[i] = dictionary_for_this_row.copy()
    df_interpolated = df_interpolated.append(df_temporary, ignore_index=True)
    plt.title(
        f"Subtrate concentrations: isocyanide {row['ic001']} M,\namine:{row['am001']} M, aldehyde:{row['ald001']} M")
    plt.xlabel('Catalyst (p-TSA) concentration (mol/L)')
    plt.ylabel('Yield (%)')
    # force the y limits to be from 0 to 0.19
    plt.ylim(0, 0.19)
    # plt.gca().yaxis.set_major_formatter(mtick.PercentFormatter())
    # vals = ax.get_yticks()
    # plt.gca().set_yticklabels(['{:}%'.format(x*100) for x in vals])
    plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda y, _: '{:.1%}'.format(y)))
    plt.subplots_adjust(left=0.15, right=0.9, top=0.85, bottom=0.15)
    plt.savefig(data_folder + run_name + f'results/interpolation_plots/interpolation_across_catalyst_{index}.png')
    if do_plot:
        plt.show()
    else:
        plt.close()

# save to file
df_interpolated.to_csv(
    data_folder + run_name + f'results/interpolated_product_concentration_together_with_model_2025-02-16.csv',
    index=False)
