from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import emcee
from _decimal import Decimal
from tqdm import tqdm
from multiprocessing import Pool
import acid_base_equilibrium as abe
import importlib
import logging

logging.basicConfig(level=logging.ERROR)

NCORES = 77
target_folder = 'visualize-results/examples/kinetics_models/'
# target_folder = ''

organize_run_results = importlib.import_module("misc-scripts.organize_run_results")


list_of_runs = tuple(['2024-03-04-run01',
                      '2024-03-04-run02'])

substances = ['c#ethyl_acetoacetate',  'c#methoxybenzaldehyde', 'c#ammonium_acetate']
substance_titles = ['Acetoacetate', 'Methoxy', 'Ammonium acetate']
# substrates = ['c#SN1OH03', 'c#HBr']

df_results = organize_run_results.join_data_from_runs([f'BPRF/{x}/' for x in list_of_runs],
                                 round_on_columns=None)
column_to_plot = 'yield#bb017'
df_results.dropna(subset=[column_to_plot], inplace=True)
df_results = df_results[~df_results[column_to_plot].isin([np.inf, -np.inf])]

# pickle df_results
df_results.to_pickle(f'{target_folder}hntz_df_results.pkl')

# unpickle df_results
# df_results = pd.read_pickle('hntz_df_results.pkl')
df_results = pd.read_pickle(f'{target_folder}hntz_df_results.pkl')

# pickle df_results
# xs = df_results['c#ethyl_acetoacetate']
# ys = df_results['c#methoxybenzaldehyde']
# plt.scatter(xs, ys)
# plt.show()

def model_of_yield_for_one_condition(*popt):
    # get row of df_results with index_in_df
    index_in_df, k1, k_1,\
        k2, k_2,\
        k3, k_3,\
        k4, k_4 = popt

    # K_eq_a = 0.0047
    # k_a = 0

    row = df_results.loc[index_in_df]
    reaction_time_in_hours = 48

    # get concentrations of substrates
    c_ald_0 = row['c#methoxybenzaldehyde']
    c_eaa_0 = row['c#ethyl_acetoacetate']
    c_aa_0 = row['c#ammonium_acetate']

    def right_hand_side_of_ode(t, y):

        # replace negative concentrations with small positive one
        concmin = 1e-12
        y_corrected = np.max((y, concmin*np.ones_like(y)), axis=0)
        c_ald, c_eaa, c_aa, c_he, c_piper, c_dm40, c_dm88 = y_corrected
        dc_ald = dc_eaa = dc_aa = dc_he = dc_piper = dc_dm40 = dc_dm88 = 0

        # 1*ALD + 2*EAA + 1*AA ⇋ HE
        forward_reaction_rate = k1 * c_ald * c_eaa**2 * c_aa
        backward_reaction_rate = k_1 * c_he
        net_reaction_rate = forward_reaction_rate - backward_reaction_rate
        dc_ald -= forward_reaction_rate
        dc_eaa -= 2 * forward_reaction_rate
        dc_aa -= forward_reaction_rate
        dc_he += net_reaction_rate

        # 2*ALD + 1*EAA + 2*AA ⇋ PIPER
        forward_reaction_rate = k2 * c_ald**2 * c_eaa * c_aa**2
        backward_reaction_rate = k_2 * c_piper
        net_reaction_rate = forward_reaction_rate - backward_reaction_rate
        dc_ald -= 2 * forward_reaction_rate
        dc_eaa -= forward_reaction_rate
        dc_aa -= 2 * forward_reaction_rate
        dc_piper += net_reaction_rate

        # 1*ALD + 1*EAA ⇋ dm40
        forward_reaction_rate = k3 * c_ald * c_eaa
        backward_reaction_rate = k_3 * c_dm40
        net_reaction_rate = forward_reaction_rate - backward_reaction_rate
        dc_ald -= forward_reaction_rate
        dc_eaa -= forward_reaction_rate
        dc_dm40 += net_reaction_rate

        # 1*ALD + 1*EAA + 1*AA ⇋ dm88_4
        forward_reaction_rate = k4 * c_ald * c_eaa * c_aa
        backward_reaction_rate = k_4 * c_dm88
        net_reaction_rate = forward_reaction_rate - backward_reaction_rate
        dc_ald -= forward_reaction_rate
        dc_eaa -= forward_reaction_rate
        dc_aa -= forward_reaction_rate
        dc_dm88 += net_reaction_rate

        res = tuple([dc_ald, dc_eaa, dc_aa, dc_he, dc_piper, dc_dm40, dc_dm88])
        res = np.array(res)

        return res



    starting_zero_concentration_filler = 1e-12
    y0 = [starting_zero_concentration_filler] * 7
    y0[0] = c_ald_0 + starting_zero_concentration_filler
    y0[1] = c_eaa_0 + starting_zero_concentration_filler
    y0[2] = c_aa_0 + starting_zero_concentration_filler
    y0 = np.array(y0)

    # USE ODE solver SOLVE_IVP to integrate from t=0 to t=14 and evaluate at t=14
    try:
        sol = solve_ivp(right_hand_side_of_ode, np.copy(np.array((0, reaction_time_in_hours))), np.copy(y0), t_eval=np.copy(tuple([reaction_time_in_hours])),
                        first_step=0.00001, max_step=8, atol=1e-8, rtol=1e-6)
        # # Old version
        # sol = solve_ivp(right_hand_side_of_ode, np.copy(np.array((0, reaction_time_in_hours))), np.copy(y0),
        #                 t_eval=np.copy(tuple([reaction_time_in_hours])),
        #                 first_step=0.00001, max_step=2, method='BDF', atol=1e-7, rtol=1e-4)

        try:
            solution = sol.y
            if type(solution) == np.ndarray:
                if solution.shape == (y0.shape[0], 1):
                    c_he_final = solution[3, 0]
                    c_piper_final = solution[4, 0]
                elif solution.shape == (y0.shape[0],):
                    c_he_final = solution[3]
                    c_piper_final = solution[4]
            elif type(solution) == list:
                if len(solution) == y0.shape[0]:
                    c_he_final = float(solution[3])
                    c_piper_final = float(solution[4])
                else:
                    print(f'Weird shape of sol.y, len={len(solution)}, type={type(solution)}. Returning zero yield.')
                    print(f'Sol status, message, success: {sol.status}, {sol.message}, {sol.success}')
                    return 0
        except TypeError:
            # print error message but not raise error
            print(f'TypeError, solution type={type(solution)}. Returning zero yield.')
            return 0

    except ValueError as e:
        # print(f'Failed to integrate at index {index_in_df}.')
        # print(f'with parameters: {popt}')
        # print(f'with y0: {y0}')
        # raise e
        c_he_final = c_piper_final = 0

    # calculate yield for HE
    least_abundant_substrate_with_stoichiometry = min(c_ald_0, c_eaa_0/2, c_aa_0)
    yield_he = c_he_final / least_abundant_substrate_with_stoichiometry

    # calculate yield for PIPER
    least_abundant_substrate_with_stoichiometry = min(c_ald_0/2, c_eaa_0, c_aa_0/2)
    yield_piper = c_piper_final / least_abundant_substrate_with_stoichiometry

    return yield_he, yield_piper


def work_wrapper(list_of_args):
    return model_of_yield_for_one_condition(*list_of_args)


def model_of_yield_for_many_conditions(*args):
    # print('Evaluating yield model for many conditions...')
    indices = args[0]
    other_params = args[1:]

    # # parallel version with pool map
    list_of_parameter_lists = [[i] + list(other_params) for i in indices]
    with Pool(NCORES) as pool:
        results = list(tqdm(pool.imap(work_wrapper, list_of_parameter_lists), total=len(indices)))
    return results


### Only the diagonal
indices_here = df_results[np.isclose(df_results['c#methoxybenzaldehyde']+df_results['c#ethyl_acetoacetate'], 110/1000)].index
# sort indices by the value of 'c#methoxybenzaldehyde'
indices_here = list(df_results.loc[indices_here].sort_values(by='c#methoxybenzaldehyde').index)

# ### All the df indices
indices_here = df_results.index

# print(f'number of indices: {len(indices_here)}')

factor = 500
p0 = np.array([factor*0.01, factor*1e-7,
               factor*0.1, factor*1e-7,
               factor*0.1/1000, factor*0.01/1000,
               factor*0.1/1000, factor*0.01/1000])

def print_parameters_with_names(theta):
    parameter_names = ['k1', 'k_1', 'k2', 'k_2', 'k3', 'k_3', 'k4', 'k_4']
    for i, name in enumerate(parameter_names):
        print(f'{name}: {theta[i]}')


# for all pkas, the lower bound is -10
lower_bounds = [0] * 8
upper_bounds = [100000] * 8

experimental_yields_he = df_results.loc[indices_here]['yield#HRP01'].to_numpy()
experimental_yields_piper = df_results.loc[indices_here]['yield#bb017'].to_numpy()


def residuals_parallel(theta):
    print(f'computing residuals in parallel for theta={theta}')
    np.save('theta_intermediate_during_opt.npy', theta)
    assert len(theta) == 8
    list_of_parameter_lists = [[i] + list(theta) for i in indices_here]
    with Pool(NCORES) as pool:
        model = list(tqdm(pool.imap(work_wrapper, list_of_parameter_lists), total=len(indices_here)))

    res = []
    for i, model_pair in enumerate(model):
        res.append(model_pair[0] - experimental_yields_he[i])
    for i, model_pair in enumerate(model):
        res.append(model_pair[1] - experimental_yields_piper[i])

    # with realistic sigmas
    sigmas = [df_results.loc[i, 'yielderr#HRP01'] for i in indices_here] + \
                [df_results.loc[i, 'yielderr#bb017'] for i in indices_here]

    sigmas = np.array(sigmas)

    res = np.array(res) / sigmas
    # sigma2 = 0.02 ** 2
    # LL = -0.5 * np.sum(res ** 2 / sigma2)
    LL = -0.5 * np.sum(res ** 2)
    print(f'Corresponding loglikelihood: {LL}')
    return res


# def plot_yields_on_one_panel(do_load=False, filename='theta_max.npy'):
#     if filename is None:
#         popt = p0
#     else:
#         popt = np.load(filename)
#
#     print(f'Optimal parameters: {popt}')
#     # plot the func
#     # local_max
#     colors = ['C0', 'C2', 'C3', 'C4', 'C5']
#     labels = ['amine min', 'amine max', 'intermediate', 'low IC low AM', 'low IC high AM']
#     if not do_load:
#         ys_all = model_of_yield_for_many_conditions(indices_here, *popt)
#     fig = plt.figure(figsize=(12, 9))
#     for i, indices in enumerate(list_of_index_sets):
#         xs = df_results.loc[indices]['ptsa']
#         if not do_load:
#             ys1 = []
#             for k, ys in enumerate(ys_all):
#                 if indices_here[k] in indices:
#                     ys1.append(ys)
#             ys2 = df_results.loc[indices]['yield']
#             np.save(f'ys1_{i}.npy', ys1)
#             np.save(f'ys2_{i}.npy', ys2)
#         else:
#             ys1 = np.load(f'ys1_{i}.npy')
#             ys2 = np.load(f'ys2_{i}.npy')
#         plt.plot(xs, ys1, '-', label=f'Model, {labels[i]}', color=colors[i])
#         plt.scatter(xs, ys2, label=f'Data, {labels[i]}', color=colors[i])
#         print(f'Max model yield for {labels[i]}: {np.max(ys1)}')
#     plt.legend()
#     plt.xlabel('Starting concentration of acid, M')
#     plt.ylabel('Yield')
#
#     fig.savefig(
#         'yieldplot.png', dpi=300)
#     # plt.show()


def plot_yields(filename='theta_max.npy', theta=None,
                do_show=False, sizefactor=5, show_legends=True, same_ylims=False):
    if not (theta is None):
        p0 = theta
    if filename is None:
        popt = p0
    else:
        popt = np.load(filename)

    print(f'Optimal parameters: {popt}')
    print_parameters_with_names(popt)


    model = model_of_yield_for_many_conditions(indices_here, *popt)
    # ys_all = []

    # largest_yield = max(np.max(ys_all), np.max(df_results.loc[indices_here]['yield']))
    # fig = plt.figure(figsize=(12, 9))
    # make number of subplots equal to the number of elements in the list_of_index_sets
    fig = plt.figure(figsize=(6, 5))
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()

    def tick_function(X):
        V = 110 - X
        return ["%.0f" % z for z in V]

    labels = ['HRP01', 'bb017']
    colors = [f'C{i}' for i in range(4)]
    colors[1] = 'gold'
    for idx in [0, 1]:
        model_ys = []
        for i, model_pair in enumerate(model):
            model_ys.append(model_pair[idx])
        model_ys = np.array(model_ys)
        experimental_ys = df_results.loc[indices_here]['yield#'+labels[idx]]
        print(f'idx {idx}, max yield: {np.max(model_ys)}')
        experimental_xs = df_results.loc[indices_here]['c#ethyl_acetoacetate']
        yerr = df_results.loc[indices_here]['yielderr#'+labels[idx]]
        ax1.errorbar(110 - experimental_xs * 1000,
                     experimental_ys * 100,
                     yerr=yerr*100, fmt='o', label=labels[idx], color=colors[idx],
                     alpha=0.5, capsize=3, capthick=1, elinewidth=1, markersize=3)
        ax1.plot(110 - experimental_xs * 1000, model_ys * 100, '-', color=colors[idx], alpha=0.5)

    ax1.set_xlim(10, 100)
    ax1Ticks = ax1.get_xticks()
    ax2Ticks = ax1Ticks
    ax2.set_xticks(ax2Ticks)
    ax2.set_xbound(ax1.get_xbound())
    ax2.set_xticklabels(tick_function(ax2Ticks))

    ax2.set_xlabel("EAA, mM")
    ax1.set_xlabel('Aldehyde, mM')

    ax1.set_ylabel('Yield, %')
    ax1.legend()
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    pass
    # theta = np.array([8.10448175e+00, 1.36644571e-03,
    #                2.45759559e+02, 1.84622923e-02,
    #                1.05885445e+00, 5.63083827e+01,
    #                5.02497025e-02, 7.90050804e+00])
    # plot_yields(filename='visualize-results/examples/kinetics_models/theta_opt.npy')

    popt = np.load(f'{target_folder}theta_opt.npy')
    # popt = np.load(f'{target_folder}theta_opt_2024-09-25.npy')
    indices_here = df_results.index
    model = model_of_yield_for_many_conditions(indices_here, *popt)
    for i, model_pair in enumerate(model):
        df_results.loc[i, 'yield#HRP01'] = model_pair[0]
        df_results.loc[i, 'yield#bb017'] = model_pair[1]

    import pickle
    pickle.HIGHEST_PROTOCOL = 4
    df_results.to_hdf('hntz_df_results_model_p4.hdf', key='df')

    # target_folder = 'visualize-results/examples/kinetics_models/'
    # params = np.load(f'{target_folder}theta_opt.npy')
    # errors = np.load(f'{target_folder}theta_perr_opt.npy')
    # param_names = ['k1', 'k_1', 'k2', 'k_2', 'k3', 'k_3', 'k4', 'k_4']
    # for i, p in enumerate(params):
    #     print(f'Kinetic rate constant {param_names[i]}: {p} +/- {errors[i]/p*100:.1f}%')

