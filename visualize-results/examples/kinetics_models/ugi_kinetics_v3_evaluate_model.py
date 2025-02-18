from ugi_kinetics_v3_emcee import *
import os

if __name__ == '__main__':
    # popt = np.load('v3_last_point_of_optimization.npy')
    # for target_id in range(5):
    #     ptsa_target = ptsa_targets[target_id]
    #     # first index where ptsa is greater than
    #     ith_ptsa = next(i for i, x in enumerate(sorted_unique_values_of_ptsa_column) if x > ptsa_target)
    #     # ith_ptsa = len(sorted_unique_values_of_ptsa_column) - 1
    #     print(f'ith_ptsa={ith_ptsa}, PTSA_value: {sorted_unique_values_of_ptsa_column[ith_ptsa]}')
    #     indices_here = df_results[df_results['ptsa'] == sorted_unique_values_of_ptsa_column[ith_ptsa]].index
    #     ys_all = model_of_yield_for_many_conditions(indices_here, *popt)
    #     # write ys_all as 'yield' column at indices_here into df_results
    #     df_results.loc[indices_here, 'yield'] = ys_all
    #
    # # save the updated df_results
    # df_results.to_csv(f'interpolated_product_concentration_model.csv', index=False)

    # make a df_results dataframe with columns 'ic001', 'am001', 'ald001', 'ptsa', 'yield'
    # Same must be uncommented in the ugi_kinetics_v3_emcee.py for this to work in parallel
    list_to_populate = []
    for ald001 in np.linspace(0.12, 0.3, 77):
        am001 = 0.3 + 0.12 - ald001
        ic001 = 0.3
        ptsa = am001
        yield_here = 0
        list_to_populate.append([ic001, am001, ald001, ptsa, yield_here])
    df_results = pd.DataFrame(list_to_populate, columns=['ic001', 'am001', 'ald001', 'ptsa', 'yield'])

    # unique_amine = list(df_results['am001'].unique())
    # sorted_amine = np.sort(unique_amine)
    # print(f'sorted_amine: {sorted_amine}')
    # yields = []
    # prods = []
    # lims = []
    # indices_to_plot = []
    # for i, amine in enumerate(sorted_amine):
    #     # if amine = 0.24
    #     df_local = df_results[np.isclose(df_results['ic001'], 0.3) & np.isclose(df_results['am001'], amine) & (
    #         np.isclose(df_results['ald001'], 0.3 + 0.12 - amine))]
    #     # find the value of ptsa closest to 'amine" value
    #     ptsa_all = df_local['ptsa'].unique()
    #     ptsa_closest_to_am = ptsa_all[np.argmin(np.abs(ptsa_all - amine))]
    #     df_local_2 = df_local[df_local['ptsa'] == ptsa_closest_to_am]
    #     indices_to_plot.append(df_local_2.index[0])
    #
    # indices_to_plot = np.array(indices_to_plot)
    # print(f'indices_to_plot: {indices_to_plot}')

    ##### YIELDS ON DIAGONAL
    popt = np.load('theta_v3_REV_opt_85c_2025-02-15.npy')
    indices_here = df_results.index
    # indices_here = indices_to_plot
    ys_all = model_of_yield_for_many_conditions(indices_here, *popt)
    df_results['yield'] = -2
    # write ys_all as 'yield' column at indices_here into df_results
    df_results.loc[indices_here, 'yield'] = ys_all

    # save the updated df_results
    df_results.to_csv(f'ugi_model_yields_on_diagonal_85c_2025-02-15.csv', index=False)

    # ###### P3/P2 ON DIAGONAL
    # popt = np.load('theta_v3_REV_opt_85c_2025-02-15.npy')
    # indices_here = df_results.index
    # # indices_here = indices_to_plot
    # ys_all = model_of_yield_for_many_conditions(indices_here, *popt)
    # # write ys_all as 'yield' column at indices_here into df_results
    # df_results.loc[indices_here, 'p3/p2'] = ys_all
    # # save the updated df_results
    # df_results.to_csv(f'ugi_model_yields_on_diagonal_p3p2.csv', index=False)

    ## ON ENTIRE RAW DATASET
    print('Evaluating on entire dataset...')
    popt = np.array([6.44617571e+03, 4.10427774e+00, 1.53526470e+04, 1.75814158e+05,
                1.07597896e+04, 3.99836093e-03, 1.93125860e+01, 3.93206588e+02,
                4.25161632e+00, 1.44080357e+03, 2.82985417e+02, 1.56419500e+02,
                9.78318054e+11, 6.17223183e+03, 4.66713746e-01, 1.24863230e+04,
                1.26205669e+04, -2.47039970e+00, 1.11688848e+01, -7.05463135e+00,
                7.28776956e+00, 8.40525819e+00, 1.11691807e+01, 1.44053613e+01,
                -3.90518238e+00, 1.66200106e+01, -5.34345656e-01, 1.45007127e+01,
                -2.50982073e-01, 1.10025081e+01, 6.70878313e+00])
    print(f'Parameters with names:')
    print_parameters_with_names(popt)
    indices_here = df_results.index
    print(f'Length of indices_here: {len(indices_here)}')
    # indices_here = indices_to_plot
    ys_all = model_of_yield_for_many_conditions(indices_here, *popt)
    # write ys_all as 'yield' column at indices_here into df_results
    df_results.loc[indices_here, 'v3_model_yield'] = ys_all

    # save the updated df_results
    df_results.to_csv(f'ugi_v3_model_yields_for_entire_exp_dataset_2025-01-31.csv', index=False)