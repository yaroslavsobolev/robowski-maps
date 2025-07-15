from robowski.settings import *
'''
This script is used to optimize the parameters of the model of the reaction. The model is a system of ODEs that describes
the kinetics of the reaction. See the Supplementary Materials Section 5.3 for the model description.
This particular script implements the model in the 'reduced' version, without the 'hydrative' mechanism.
'''

from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import emcee
from _decimal import Decimal
from tqdm import tqdm
from multiprocessing import Pool
import robowski.kinetics_models.acid_base_equilibrium as abe


# experiment_name = 'multicomp-reactions/2023-06-19-run01/'
df_results = pd.read_csv('interpolated_product_concentration.csv')
df_results['yield'] = df_results['yield'].apply(lambda x: 0 if x < 0 else x)
substances = ['ic001','am001','ald001','ptsa']
substance_titles = ['Isocyanide', 'Amine', 'Aldehyde', 'p-TSA']

padding_rows_count = (df_results[substances] == 0).all(axis=1).sum()
# print(f"There are {padding_rows_count} padding rows (with zero concentrations of substrates).")
df_results = df_results[(df_results[substances] != 0).any(axis=1)]
substrate_cs = []
for substance in substances:
    substrate_cs.append(df_results[substance].to_numpy())
xs0, ys0, zs0, cats = substrate_cs
# print('Max concentrations of substrates: ')
# for x in [xs0, ys0, zs0]:
#     print(max(x))
minimal_concentration_of_substrates = np.min(np.array([xs0, ys0, zs0]))
maximum_concentration_of_substrates = np.max(np.array([xs0, ys0, zs0]))
unique_cats = list(sorted(list(df_results['ptsa'].unique())))

def model_of_yield_for_one_condition(*popt):
    # get row of df_results with index_in_df

    theta_0 = np.array([  3.29816543e+02, 3.47572969e+01,  # k1, k_1, \
                        3.13158495e+03, 100*7.72938177e+02,  # k2, k_2, \
                        100*5.08931575e+02, 3.90934152e-03,  # k3, k_3, \
                        100*2.05835089e-01, 8.03029698e+00,  # k4, k_4, \
                        8.74858657e+01, # k5, \
                        0, 0, # k6, k_6, \
                        0,      # k7, \
                        0,      # k8, \
                        0,# k9, \
                        0.2, 0,      # km1, k_m1, \
                        0,  # km3, \
                        0.004499,# K_eq_a, \
                        0, # k_a, \
                        -1.7, # pKa_ptsa, \ # -1.9257310211345142
                        12.23, # pKa_BuNH3, \ # 12.9705271
                        -7.2, # pKa_DMFH, \
                        8, # pKa_p1, \ # 7.7645547494439295
                        6.1, # pKa_p2, \
                        10, # pKa_p3, \
                        13.1, # pKa_c, \
                        0.86, # pKa_ch, \
                        17, # pKa_q1, \
                        0, # pKa_q2h, \
                        17, # pKa_m1h, \
                        0, # pKa_m4h, \
                        11.97, # pKa_me2nh2, \
                        5.9 # pKa_q4
                      ])

    k1, k_1,\
        k2, k_2,\
        k3, k_3,\
        k4, k_4,\
        k5,\
        k6, k_6,\
        k7,\
        k8,\
        k9,\
        km1, k_m1,\
        km3,\
        K_eq_a,\
        k_a,\
        pKa_ptsa,\
        pKa_BuNH3,\
        pKa_DMFH,\
        pKa_p1,\
        pKa_p2,\
        pKa_p3,\
        pKa_c,\
        pKa_ch,\
        pKa_q1,\
        pKa_q2h,\
        pKa_m1h,\
        pKa_m4h,\
        pKa_me2nh2,\
        pKa_q4 = theta_0

    index_in_df, k1, k_1,\
        k2, k_2,\
        k3, k_3,\
        k4, k_4,\
        k5,\
        km1, k_m1,\
        pKa_ptsa,\
        pKa_BuNH3,\
        pKa_p1 = popt



    row = df_results.loc[index_in_df]
    npoints = 100
    reaction_time_in_hours = 14 # when fitting experimental data obtained at 16 hours reaction time, we leave this variable equal to 14 hours and then scale the kinetic rate constants by 14/16 when reporting them in the article. Mathematically, this scaling is equivalent to setting this parameter to 16 hours here in the code; we use clumsy scaling approach for legacy reasons.
    dt = reaction_time_in_hours / npoints
    # get concentrations of substrates
    c_a = row['ald001']
    c_b = row['am001']
    c_c = row['ic001']
    c_ptsa = row['ptsa']
    initial_concentration_of_least_abundant_susbtrate = min(c_c, c_b, c_a)

    # initial concentrations
    # c_p1 = 0
    # c_p2 = 0
    # c_salt = 0
    # c_product = 0
    # c_h2o = c_ptsa
    c_dmf = 12.915

    def right_hand_side_of_ode(t, y):

        # replace negative concentrations with small positive one
        concmin = 1e-12
        y_corrected = np.max((y, concmin*np.ones_like(y)), axis=0)
        c_a, c_b, c_c, c_p1, c_p2, c_p3, c_p4, c_me2nh2, c_q1, c_q2, c_q4, c_m1, c_m4, c_h2o, c_a_h2o, c_product = y_corrected
        dc_a = dc_b = dc_c = dc_p1 = dc_p2 = dc_p3 = dc_p4 = dc_me2nh2 = dc_q1 = dc_q2 = dc_q4 = dc_m1 = dc_m4 = dc_h2o = dc_a_h2o = dc_product = 0

        # Figuring our the pH and the concentrations of different ions
        substances = ({'pKa': pKa_BuNH3, 'conc': c_b, 'charge': 1},  # ammonium
                      {'pKa': pKa_ptsa, 'conc': c_ptsa, 'charge': 0}, # TsOH + H2O -> TSO- + H3O+
                      {'pKa': pKa_DMFH, 'conc': c_dmf, 'charge': 1}, # DHFH+ + H2O -> DMF + H3O+
                      {'pKa': pKa_p1, 'conc': c_p1, 'charge': 1},
                      {'pKa': pKa_p2, 'conc': c_p2, 'charge': 1},
                      {'pKa': pKa_p3, 'conc': c_p3, 'charge': 1},
                      {'pKa': pKa_ch, 'pKa2': pKa_c, 'conc': c_c, 'charge': 1},
                      {'pKa': pKa_q1, 'conc': c_q1, 'charge': 0},
                      {'pKa': pKa_q2h, 'conc': c_q2, 'charge': 1},
                      {'pKa': pKa_m1h, 'conc': c_m1, 'charge': 1},
                      {'pKa': pKa_m4h, 'conc': c_m4, 'charge': 1},
                      {'pKa': pKa_me2nh2, 'conc': c_me2nh2, 'charge': 1},
                      {'pKa': pKa_q4, 'conc': c_q4, 'charge': 1}
                      )
        # convert all 'pKa' and 'conc' to Decimal for speed
        # for sub in substances:
        #     sub['pKa'] = Decimal(sub['pKa'])
        #     sub['conc'] = Decimal(sub['conc'])
        # Kw_here = 10 ** (-27.44) * (13.9) ** 2

        # Entry 4 in Table 1 in the article [Clara Rifols, 1994, "Autoprotolysis in aqueous organic solvent mixtures
        # Water-amide and water-amine binary systems"] gives autoprotolysis of H2O in DMF
        Kw_here = Decimal(10) ** Decimal(-27.44) / (abe.PURE_WATER_MOLARITY) ** 2

        try:
            ph_here = abe.solve_for_zero_charge(water_concentration=c_h2o, kW=Kw_here, substances=substances)
        except Exception as e:
            print(f'Failed to solve for zero charge at index {index_in_df}.')
            print(f'with parameters: {popt}')
            print(f'with y: {y_corrected}')
            raise e
        subs, remaining_water, oh_minus = abe.concentrations_by_ph(water_concentration=c_h2o, Kw=Kw_here,
                                                               substances=substances,
                                                               ph=ph_here, return_Decimals=False)
        eq_BuNH3 = subs[0]
        eq_ptsa = subs[1]
        eq_DMFH = subs[2]
        eq_p1 = subs[3]
        eq_p2 = subs[4]
        eq_p3 = subs[5]
        eq_ch = subs[6]
        eq_q1 = subs[7]
        eq_q2h = subs[8]
        eq_m1h = subs[9]
        eq_m4h = subs[10]
        eq_me2nh2 = subs[11]
        eq_q4 = subs[12]

        # A+B ⇋ P_1
        forward_reaction_rate = k1 * c_a * eq_BuNH3['conc_deprot']
        backward_reaction_rate = k_1 * eq_p1['conc_deprot']
        net_reaction_rate = forward_reaction_rate - backward_reaction_rate
        dc_p1 += net_reaction_rate
        dc_a += -net_reaction_rate
        dc_b += -net_reaction_rate

        # P_1^+ ⇋ P_2^+ + H2O
        forward_reaction_rate = k2 * eq_p1['conc_prot']
        backward_reaction_rate = k_2 * eq_p2['conc_prot'] * remaining_water
        net_reaction_rate = forward_reaction_rate - backward_reaction_rate
        dc_p1 += -net_reaction_rate
        dc_p2 += net_reaction_rate
        dc_h2o += net_reaction_rate

        # P_2^+ + B ⇋ P_3^+
        forward_reaction_rate = k3 * eq_p2['conc_prot'] * eq_BuNH3['conc_deprot']
        backward_reaction_rate = k_3 * eq_p3['conc_prot']
        net_reaction_rate = forward_reaction_rate - backward_reaction_rate
        dc_p2 += -net_reaction_rate
        dc_b += -net_reaction_rate
        dc_p3 += net_reaction_rate

        # P_2^+ + C ⇋ P_4+
        forward_reaction_rate = k4 * eq_p2['conc_prot'] * eq_ch['conc_deprot']
        backward_reaction_rate = k_4 * c_p4
        net_reaction_rate = forward_reaction_rate - backward_reaction_rate
        dc_p2 += -net_reaction_rate
        dc_c += -net_reaction_rate
        dc_p4 += net_reaction_rate

        # P_4+DMF ⇋ Prod + Me2NH2
        forward_reaction_rate = k5 * c_p4
        backward_reaction_rate = 0
        net_reaction_rate = forward_reaction_rate - backward_reaction_rate
        dc_p4 += -net_reaction_rate
        dc_product += net_reaction_rate
        dc_me2nh2 += net_reaction_rate

        # C_minus + A ⇋ Q_1minus
        forward_reaction_rate = k6 * eq_ch['conc_deprot_2'] * c_a
        backward_reaction_rate = k_6 * eq_q1['conc_deprot']
        net_reaction_rate = forward_reaction_rate - backward_reaction_rate
        dc_c += -net_reaction_rate
        dc_a += -net_reaction_rate
        dc_q1 += net_reaction_rate

        # Q_1 → Q_4
        forward_reaction_rate = k7 * eq_q1['conc_prot']
        backward_reaction_rate = 0
        net_reaction_rate = forward_reaction_rate - backward_reaction_rate
        dc_q1 += -net_reaction_rate
        dc_q4 += net_reaction_rate

        # CH+ + H2O → Q2H+
        forward_reaction_rate = k8 * eq_ch['conc_prot'] * remaining_water
        backward_reaction_rate = 0
        net_reaction_rate = forward_reaction_rate - backward_reaction_rate
        dc_c += -net_reaction_rate
        dc_h2o += -net_reaction_rate
        dc_q2 += net_reaction_rate

        # Q_4+ + H2O → Q_5
        forward_reaction_rate = k9 * eq_q4['conc_prot'] * remaining_water
        net_reaction_rate = forward_reaction_rate
        dc_q4 += -net_reaction_rate
        dc_h2o += -net_reaction_rate

        # C+A ⇋ M1
        forward_reaction_rate = km1 * eq_ch['conc_deprot'] * c_a
        backward_reaction_rate = k_m1 * eq_m1h['conc_deprot']
        net_reaction_rate = forward_reaction_rate - backward_reaction_rate
        dc_c += -net_reaction_rate
        dc_a += -net_reaction_rate
        dc_m1 += net_reaction_rate

        # M1H^+ + H2O ⇋ M4H^+
        forward_reaction_rate = km3 * eq_m1h['conc_prot'] * remaining_water
        net_reaction_rate = forward_reaction_rate
        dc_m1 += -net_reaction_rate
        dc_h2o += -net_reaction_rate
        dc_m4 += net_reaction_rate
        
        # A + H2O ⇋ A_H2O
        k_a_h2o_forward = k_a * K_eq_a
        k_a_h2o_backward = k_a
        forward_reaction_rate = k_a_h2o_forward * c_a * remaining_water
        backward_reaction_rate = k_a_h2o_backward * c_a_h2o
        net_reaction_rate = forward_reaction_rate - backward_reaction_rate
        dc_a += -net_reaction_rate
        dc_h2o += -net_reaction_rate
        dc_a_h2o += net_reaction_rate

        res = tuple([dc_a, dc_b, dc_c, dc_p1, dc_p2, dc_p3, dc_p4, dc_me2nh2, dc_q1, dc_q2, dc_q4, dc_m1, dc_m4, dc_h2o, dc_a_h2o, dc_product])
        res = np.array(res)

        return res



    starting_zero_concentration_filler = 1e-12
    y0 = [starting_zero_concentration_filler] * 16
    y0[0] = c_a + starting_zero_concentration_filler
    y0[1] = c_b + starting_zero_concentration_filler
    y0[2] = c_c + starting_zero_concentration_filler
    y0[14] = c_ptsa + starting_zero_concentration_filler  # because initial water concentration is the same as concentration of p-TSA, by definition of experiment
    y0 = np.array(y0)

    # USE ODE solver SOLVE_IVP to integrate from t=0 to t=14 and evaluate at t=14
    try:
        # sol = solve_ivp(right_hand_side_of_ode, np.copy(np.array((0, reaction_time_in_hours))), np.copy(y0), t_eval=np.copy(tuple([reaction_time_in_hours])),
        #                 first_step=0.00001, max_step=2, atol=1e-8, rtol=1e-6)

        # LSODA VERSION
        sol = solve_ivp(right_hand_side_of_ode, np.copy(np.array((0, reaction_time_in_hours))), np.copy(y0), t_eval=np.copy(tuple([reaction_time_in_hours])),
                        first_step=0.00001, max_step=2, atol=1e-8, rtol=1e-6, method='LSODA')
        # # Old version
        # sol = solve_ivp(right_hand_side_of_ode, np.copy(np.array((0, reaction_time_in_hours))), np.copy(y0),
        #                 t_eval=np.copy(tuple([reaction_time_in_hours])),
        #                 first_step=0.00001, max_step=2, method='BDF', atol=1e-7, rtol=1e-4)

        try:
            solution = sol.y
            if type(solution) == np.ndarray:
                if solution.shape == (y0.shape[0], 1):
                    c_product = solution[-1, 0]
                elif solution.shape == (y0.shape[0],):
                    c_product = solution[-1]
            elif type(solution) == list:
                if len(solution) == y0.shape[0]:
                    c_product = float(solution[-1])
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
        c_product = 0

    # def right_hand_side_of_ode_z(t, z):
    #     # assert that np.exp(z) is positive and finite
    #     if np.any(np.isnan(np.exp(z))) or np.any(np.isinf(np.exp(z))):
    #         print(f'exp(z) is nan or inf: {np.exp(z)}')
    #         raise ValueError('z is nan or inf')
    #     return np.array(right_hand_side_of_ode(t, np.exp(z))) / np.exp(z)
    # # here I tried to do it with logarithmic derivatives -- did not work
    # # USE ODE solver SOLVE_IVP to integrate from t=0 to t=14 and evaluate at t=14
    # sol = solve_ivp(right_hand_side_of_ode_z, [0, reaction_time_in_hours], np.log(y0), t_eval=[reaction_time_in_hours],
    #                 first_step=0.00000001, max_step=0.1)
    # c_product = np.exp(sol.y[-1, 0])

    reaction_yield = c_product / initial_concentration_of_least_abundant_susbtrate
    # print(f'Evaluated yield for index {index_in_df}: {reaction_yield}')
    # print(f'for parameters: {popt}')
    return reaction_yield

def work_wrapper(list_of_args):
    return model_of_yield_for_one_condition(*list_of_args)

def model_of_yield_for_many_conditions(*args):
    # print('Evaluating yield model for many conditions...')
    indices = args[0]
    other_params = args[1:]
    # # serial version
    # return [model_of_yield_for_one_condition(*([i] + list(other_params))) for i in tqdm(indices)]
    # # parallel version with pool map
    list_of_parameter_lists = [[i] + list(other_params) for i in indices]
    with Pool(80) as pool:
        results = list(tqdm(pool.imap(work_wrapper, list_of_parameter_lists), total=len(indices)))
    return results


# decimate_by = 3
decimate_by = 6

list_of_index_sets = [df_results[(df_results['ic001'] == 0.3) & (df_results['am001'] == 0.12) & (df_results['ald001'] == 0.3)].index[::decimate_by],
                      df_results[(df_results['ic001'] == 0.3) & (df_results['am001'] == 0.3) & (df_results['ald001'] == 0.12)].index[::decimate_by],
                      df_results[(df_results['ic001'] == 0.3) & (df_results['am001'] == 0.21) & (df_results['ald001'] == 0.21)].index[::decimate_by],
                      df_results[(df_results['ic001'] == 0.12) & (df_results['am001'] == 0.12) & (df_results['ald001'] == 0.3)].index[::decimate_by],
                      df_results[(df_results['ic001'] == 0.12) & (df_results['am001'] == 0.3) & (df_results['ald001'] == 0.12)].index[::decimate_by]
                      ]

# concatenate indices
ptsa_ranges = [[0.04, 0.0483],
               [0.12, 0.124],
               [0.1718, 0.1739],
               [0.2368, 0.2455],
               [0.2983, 0.3]]

# # # only use indices where ptsa is in ptsa ranges
# for i, indices_list in enumerate(list_of_index_sets):
#     new_index = []
#     for idx in indices_list:
#         for ptsa_range in ptsa_ranges:
#             if ptsa_range[0] <= df_results.loc[idx]['ptsa'] <= ptsa_range[1]:
#                 new_index.append(idx)
#                 break
#     list_of_index_sets[i] = new_index[:]

indices_here = np.concatenate(list_of_index_sets)

print(f'number of indices: {len(indices_here)}')

p0 = np.array([  3.29816543e+02, 3.47572969e+01,  # k1, k_1, \
                        3.13158495e+03, 100*7.72938177e+02,  # k2, k_2, \
                        100*5.08931575e+02, 3.90934152e-03,  # k3, k_3, \
                        100*2.05835089e-01, 8.03029698e+00,  # k4, k_4, \
                        8.74858657e+01, # k5, \
                        6.67323340e-01, 6.67417607e+02,      # km1, k_m1, \
                        -1.09674313e+00, # pKa_ptsa, \
                        1.10722293e+01, # pKa_BuNH3, \ # 12.842106004719119
                        8.27868680e+00, # pKa_p1, \
                      ])

def print_parameters_with_names(theta):
    parameter_names = ['k1', 'k_1', 'k2', 'k_2', 'k3', 'k_3', 'k4', 'k_4', 'k5', 'km1', 'k_m1', 'pKa_ptsa', 'pKa_BuNH3', 'pKa_p1']
    for i, name in enumerate(parameter_names):
        print(f'{name}: {theta[i]}')


# for all pkas, the lower bound is -10
lower_bounds = [0] * 11    + [-8] * 3
upper_bounds = [5000] * 11 + [20] * 3
lower_bounds[11] = -2.5 # pKa_ptsa,\
upper_bounds[11] = 0 # pKa_ptsa,\
lower_bounds[12] = 9 # pKa_BuNH3,\
upper_bounds[12] = 15.5 # pKa_BuNH3,\
lower_bounds[13] = -1 # pKa_p1,\
upper_bounds[13] = 12 # pKa_p1,\

upper_bounds[3] = 150000 # k_2
upper_bounds[4] = 150000 # k3

experimental_yields = df_results.loc[indices_here]['yield']
ptsa_values = df_results.loc[indices_here]['ptsa']


# # plt.plot(ptsa_values, model, 'r-', label='Model')
# plt.scatter(ptsa_values, experimental_yields, label='Data')
# plt.show()

def log_likelihood(theta):
    assert len(theta) == 33
    model = [model_of_yield_for_one_condition(*([i] + list(theta))) for i in indices_here]
    sigma2 = 0.02 ** 2
    return -0.5 * np.sum((experimental_yields - model) ** 2 / sigma2)


def log_likelihood_parallel(theta):
    print(f'computing log likelihood in parallel for theta={theta}')
    np.save('theta_intermediate_during_v3bREV_opt.npy', theta)
    assert len(theta) == 14
    list_of_parameter_lists = [[i] + list(theta) for i in indices_here]
    with Pool(77) as pool:
        model = list(tqdm(pool.imap(work_wrapper, list_of_parameter_lists), total=len(indices_here)))
    sigma2 = 0.02 ** 2
    res = -0.5 * np.sum((experimental_yields - model) ** 2 / sigma2)
    print(f'Loglikelihood: {res}')
    return res

def residuals_parallel(theta):
    print(f'computing residuals in parallel for theta={theta}')
    np.save('theta_intermediate_during_v3bREV_opt.npy', theta)
    assert len(theta) == 14
    list_of_parameter_lists = [[i] + list(theta) for i in indices_here]
    with Pool(77) as pool:
        model = list(tqdm(pool.imap(work_wrapper, list_of_parameter_lists), total=len(indices_here)))
    res = experimental_yields - model

    sigma2 = 0.02 ** 2
    LL = -0.5 * np.sum((experimental_yields - model) ** 2 / sigma2)
    print(f'Corresponding loglikelihood: {LL}')
    return res

def log_prior(theta):
    if all([lower_bounds[i] < theta[i] < upper_bounds[i] for i in range(len(theta))]):
        return 0
    return -np.inf

def log_probability(theta):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta)


def plot_yields_on_one_panel(do_load=False, filename='theta_max.npy'):
    if filename is None:
        popt = p0
    else:
        popt = np.load(filename)

    print(f'Optimal parameters: {popt}')
    # plot the func
    # local_max
    colors = ['C0', 'C2', 'C3', 'C4', 'C5']
    labels = ['amine min', 'amine max', 'intermediate', 'low IC low AM', 'low IC high AM']
    if not do_load:
        ys_all = model_of_yield_for_many_conditions(indices_here, *popt)
    fig = plt.figure(figsize=(12, 9))
    for i, indices in enumerate(list_of_index_sets):
        xs = df_results.loc[indices]['ptsa']
        if not do_load:
            ys1 = []
            for k, ys in enumerate(ys_all):
                if indices_here[k] in indices:
                    ys1.append(ys)
            ys2 = df_results.loc[indices]['yield']
            np.save(f'ys1_{i}.npy', ys1)
            np.save(f'ys2_{i}.npy', ys2)
        else:
            ys1 = np.load(f'ys1_{i}.npy')
            ys2 = np.load(f'ys2_{i}.npy')
        plt.plot(xs, ys1, '-', label=f'Model, {labels[i]}', color=colors[i])
        plt.scatter(xs, ys2, label=f'Data, {labels[i]}', color=colors[i])
        print(f'Max model yield for {labels[i]}: {np.max(ys1)}')
    plt.legend()
    plt.xlabel('Starting concentration of acid, M')
    plt.ylabel('Yield')
    
    fig.savefig(
        'yieldplot.png', dpi=300)
    # plt.show()


def plot_yields(do_load=False, filename='theta_max.npy', theta=None, do_show=False, sizefactor=5, show_legends=True, same_ylims=False):
    if not (theta is None):
        p0 = theta
    if filename is None:
        popt = p0
    else:
        popt = np.load(filename)

    print(f'Optimal parameters: {popt}')
    print_parameters_with_names(popt)
    # plot the func
    # local_max
    colors = ['C0', 'C2', 'C3', 'C4', 'C5']
    labels = ['amine min', 'amine max', 'intermediate', 'low IC low AM', 'low IC high AM']
    if not do_load:
        ys_all = model_of_yield_for_many_conditions(indices_here, *popt)
        print(f'Smallest yield value: {np.min(ys_all)}')

    largest_yield = max(np.max(ys_all), np.max(df_results.loc[indices_here]['yield']))
    # fig = plt.figure(figsize=(12, 9))
    # make number of subplots equal to the number of elements in the list_of_index_sets
    fig, axarr = plt.subplots(len(list_of_index_sets) // 2 + 1, 2, figsize=(2*sizefactor, (len(list_of_index_sets) // 2 + 1)*sizefactor), sharex=True)
    # flatten the axarr list
    axarr[-1, 0].set_xlabel('Starting concentration of acid, M')
    axarr[-1, 1].set_xlabel('Starting concentration of acid, M')
    axarr = [item for sublist in axarr for item in sublist]
    for i, indices in enumerate(list_of_index_sets):
        xs = df_results.loc[indices]['ptsa']
        if not do_load:
            ys1 = []
            for k, ys in enumerate(ys_all):
                if indices_here[k] in indices:
                    ys1.append(ys)
            ys2 = df_results.loc[indices]['yield']
            np.save(f'ys1_{i}.npy', ys1)
            np.save(f'ys2_{i}.npy', ys2)
        else:
            ys1 = np.load(f'ys1_{i}.npy')
            ys2 = np.load(f'ys2_{i}.npy')
        axarr[i].plot(xs, ys1, '-', label=f'Model, {labels[i]}', color=colors[i])
        axarr[i].scatter(xs, ys2, label=f'Data, {labels[i]}', color=colors[i])
        if same_ylims:
            axarr[i].set_ylim([0, largest_yield])

        print(f'Max model yield for {labels[i]}: {np.max(ys1)}')
        # add legend in axarr[i]
        if show_legends:
            axarr[i].legend()
        axarr[i].set_ylabel('Yield')

    # plt.legend()
    # axarr[].set_xlabel('Starting concentration of acid, M')
    # plt.ylabel('Yield')
    plt.tight_layout()
    fig.savefig(
        'yieldplot.png', dpi=300)
    if do_show:
        plt.show()

if __name__ == '__main__':
    # plot_yields()

    # set random seed
    # np.random.seed(42)

    # # server settings
    filename = "emcee_backend_gauss_2024_08_16a.h5"
    mode = 'continue'
    nwalkers_here = 80*6
    nsteps = 200 # 10 steps take about 1 hour
    ncores = 80

    # # local settings
    # mode = 'start'
    # nwalkers_here = 70
    # nsteps = 3
    # ncores = 4

    pos = np.array(p0) * (1 + 0.9 * np.random.randn(nwalkers_here, len(p0)))
    pos[pos<0] = 0.001
    # for pKas, sample uniformly from p0-0.2 to p0+0.2
    for i in range(11, 14):
        pos[:, i] = np.random.uniform(low=p0[i] - 0.6, high=p0[i] + 0.6, size=nwalkers_here)

    # sample uniformy from ranges set by lower bounds and upper bounds
    # pos = np.array([np.random.uniform(low=lower_bounds, high=upper_bounds) for i in range(nwalkers_here)])
    nwalkers, ndim = pos.shape
    print(f'Number of walkers: {nwalkers}\n'
          f'Number of dimensions: {ndim}\n'
          f'Number of steps: {nsteps}'      )

    # Set up the backend.
    # Don't forget to clear it in case the file already exists.
    backend  = emcee.backends.HDFBackend(filename)
    if mode == 'start':
        backend.reset(nwalkers, ndim)
    print("Initial size of backend: {0}".format(backend.iteration))


    with Pool(ncores) as pool:
        sampler = emcee.EnsembleSampler(
            nwalkers, ndim, log_probability, pool=pool, backend=backend
            # moves=[
            #     (emcee.moves.RedBlueMove(),   0.4),
            #     (emcee.moves.DEMove(),        0.4),
            #     (emcee.moves.DESnookerMove(), 0.2),
            # ],
        )
        if mode == 'continue':
            sampler.run_mcmc(None, nsteps, progress=True, skip_initial_state_check=True, store=True)
        elif mode == 'start':
            sampler.run_mcmc(pos, nsteps, progress=True, skip_initial_state_check=True, store=True)

    print("Final size of backend: {0}".format(backend.iteration))

    # # nplots = 5
    # # fig, axes = plt.subplots(nplots, figsize=(10, 7), sharex=True)
    # samples = sampler.get_chain()
    # # # labels = range(len(p0))
    # # for i in range(nplots):
    # #     ax = axes[i]
    # #     ax.plot(samples[:, :, i], "k", alpha=0.3)
    # #     ax.set_xlim(0, len(samples))
    # #     # ax.set_ylabel(labels[i])
    # #     ax.yaxis.set_label_coords(-0.1, 0.5)
    # #
    # # axes[-1].set_xlabel("step number")
    # # plt.show()
    #
    # tau = sampler.get_autocorr_time()
    # print(f'tau: {tau}')
    #
    # flat_samples = sampler.get_chain(discard=30, thin=5, flat=True)
    # print(flat_samples.shape)
    #
    # best_index = np.argmax(sampler.flatlnprobability)
    # theta_max = sampler.flatchain[best_index]
    # np.save('theta_best.npy', theta_max)
    #
    # print(f'Best theta: {theta_max}')
    #
    #
    # import corner
    #
    # labels = [str(i) for i in range(len(p0))]
    # fig = corner.corner(
    #     flat_samples, labels=labels, truths=theta_max
    # )
    #
    # fig.savefig('tempcorner.png', dpi=300)

