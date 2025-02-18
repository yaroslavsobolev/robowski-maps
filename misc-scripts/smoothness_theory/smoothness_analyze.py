from smoothness_v2 import *

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

def analyze_one_reaction(interesting_reaction, target_folder = f'{data_folder}simulations/smoothness'):
    y0 = np.zeros(10)
    re2 = pickle.load(open(f'{target_folder}/{interesting_reaction}/reaction_network.pkl', 'rb'))

    ts = np.logspace(-3, 3.5, 200)
    solution = re2.solve_ode(ts, y0, first_step=0.000000011, max_step=0.1, atol=1e-10, rtol=1e-9)
    concentrations = solution.y.T
    times = solution.t

    # sum of concentration at the end of the simulation
    sum_concentrations = np.sum(concentrations[-1])
    print(f'Sum of concentrations at the end of the simulation: {sum_concentrations}')

    plt.plot(solution.t, solution.y.T, label=range(re2.N))

    # make log scale in x axis
    plt.xscale('log')
    plt.xlabel('Time')
    plt.ylabel('Concentration')
    plt.legend()
    plt.show()

    # show k-matrix
    kmat = re2.k_tensor

    # evaluating partial derivatives of elements of y at the end of the simulation by the initial concentrations y0
    partial_derivs = []
    delta_relative = 0.0001
    for i in range(re2.N):
        y0_plus_delta = y0.copy()
        delta = y0_plus_delta[i] * delta_relative
        y0_plus_delta[i] += delta
        solution_plus_delta = re2.solve_ode(ts, y0_plus_delta, first_step=0.000000011, max_step=0.1, atol=1e-10, rtol=1e-9)
        y_deriv = (solution_plus_delta.y.T - solution.y.T) / delta
        partial_derivs.append(y_deriv.copy())
        # plt.plot(solution.t, y_deriv, label=i)

    partial_derivs = np.array(partial_derivs)

    # plot the partial derivative against time
    for i in range(re2.N):
        for j in range(re2.N):
            if j == 9:
                plt.plot(solution.t, partial_derivs[i, :, j], label=f'{i} {j}')
    plt.legend()
    plt.show()
    # highest partial derivative for each given time
    max_partial_derivs = np.max(np.abs(partial_derivs), axis=(0, 2))
    print(f'Maximum partial derivative: {np.max(max_partial_derivs)}')


def make_hists_for_C_cases(case ='C1'):
    target_folder = f'{data_folder}simulations/smoothness'
    # short_db = load_short_db(suffix='_C4')
    # short_db = load_short_db(suffix='_C1_C2')
    short_db = load_short_db(suffix='_C12345')
    # short_db = load_short_db(suffix='_v2.pkl')

    # ts = np.logspace(-3, 3.5, 200)
    # for key, value in short_db.items():
    #     plt.plot(ts, value['max_over_t'], color='k', alpha=0.05)
    #
    # plt.xscale('log')
    # plt.xlabel('Time')
    # plt.ylabel('Maximum $D_(ij)$ over all initial conditions')
    # plt.show()

    # plot the max_over_t for all reactions that are weakly ergodic, strict backward and single use
    ts = np.logspace(-3, 3.5, 200)
    max_over_ts = []
    for key, value in short_db.items():
        if not (value['case'] == case):
            continue
        # if value['strict_backward'] and (value['single_use']) and (value['max_simple_cycle_length'] <= 2) and (
        # not value['weakly_ergodic_actually']):
        max_over_ts.append(value['max_over_t'])
        if np.max(value['max_over_t']) > 7:
            print(key)
    print(f'Found {len(max_over_ts)} reactions that match the criteria')
    max_over_ts = np.array(max_over_ts)
    max_over_ts_global = np.max(max_over_ts, axis=0)
    max_for_each_reaction = np.max(max_over_ts, axis=1)

    print(f'Max of max_over_ts_global: {np.max(max_over_ts_global):.8f}')
    print(f'Median of max_for_each_reaction: {np.median(max_for_each_reaction):.8f}')

    # print the top 10 reactions with the highest max_for_each_reaction. Sort then print first 10
    indices = np.argsort(max_for_each_reaction)
    print('Top 10 reactions with the highest max_for_each_reaction:')
    for i in range(10):
        print(f'{max_for_each_reaction[indices[-i-1]]}')
    print('Top 10 with commas:')
    s = ''
    for i in range(10):
        s += f'{max_for_each_reaction[indices[-i-1]]:.4f}, '
    print(s)

    f1 = plt.figure(figsize=(5, 1.5))
    plt.hist(max_for_each_reaction, bins=100, density=True, color='grey')
    plt.axvline(np.max(max_over_ts_global), color='red')
    plt.xlabel('Highest partial derivative for a given reaction')
    plt.ylabel('Probability\ndensity')
    plt.tight_layout()
    f1.savefig(f'misc-scripts/smoothness_theory/figures/hist_{case}.png', dpi=300)
    plt.show()

    # interesting_reaction = 'reaction4ffc304e-0d60-4332-83d5-93f7f2bd0377'
    # re2 = pickle.load(open(f'{target_folder}/{interesting_reaction}/reaction_network.pkl', 'rb'))
    # for condition_id, condition_scan in enumerate(re2.condition_scans):
    #     if np.max(condition_scan['max_partial_derivs_over_t']) > threshold:
    #         print(f'Found condition with max partial derivative over t > {threshold}, the index is {condition_id}')
    #         y0 = condition_scan['y0']
    # print(y0)

def plot_hist(short_db_filename, filter_callable, output_filename):
    with open(short_db_filename, 'rb') as f:
        short_db = pickle.load(f)

    # plot the max_over_t for all reactions that are weakly ergodic, strict backward and single use
    ts = np.logspace(-3, 3.5, 200)
    max_over_ts = []
    for key, value in short_db.items():
        if filter_callable(value):
            max_over_ts.append(value['max_over_t'])
            # if np.max(value['max_over_t']) > 7.83:
            #     print(key)
    print(f'Found {len(max_over_ts)} reactions that match the criteria')
    max_over_ts = np.array(max_over_ts)
    max_over_ts_global = np.max(max_over_ts, axis=0)
    max_for_each_reaction = np.max(max_over_ts, axis=1)
    print(f'Max of max_over_ts_global: {np.max(max_over_ts_global):.8f}')
    print(f'Median of max_for_each_reaction: {np.median(max_for_each_reaction):.3f}')
    print(f'{np.max(max_over_ts_global):.3f}\t{np.median(max_for_each_reaction):.3f}')

    f1 = plt.figure(figsize=(5, 1.8))
    # set right margin to zero
    plt.subplots_adjust(right=0.99)
    # set bottom margin to larger value
    plt.subplots_adjust(bottom=0.25)
    plt.hist(max_for_each_reaction, bins=100, density=True, color='grey')
    plt.axvline(np.max(max_over_ts_global), color='red')
    plt.xlabel('Highest partial derivative')
    plt.ylabel('Probability density')
    # plt.tight_layout()
    f1.savefig(f'misc-scripts/smoothness_theory/figures/{output_filename}.png', dpi=300)
    plt.show()
    # print(f'Maximum value of max_over_ts_global: {np.max(max_over_ts_global)}')

def make_hists_for_A_cases():
    plot_hist(short_db_filename=f'{data_folder}simulations/smoothness/short_db_noback.pkl',
              filter_callable=lambda value: value['single_use'],
              output_filename='hist_A1')

    plot_hist(short_db_filename=f'{data_folder}simulations/smoothness/short_db_noback.pkl',
              filter_callable=lambda value: (not value['single_use'] and value['weakly_ergodic_actually']),
              output_filename='hist_A2')

    plot_hist(short_db_filename=f'{data_folder}simulations/smoothness/short_db_noback.pkl',
              filter_callable=lambda value: (not value['single_use'] and (not value['weakly_ergodic_actually'])),
              output_filename='hist_A3')

    plot_hist(short_db_filename=f'{data_folder}simulations/smoothness/short_db_manycycles.pkl',
              filter_callable=lambda value: (value['single_use'] and (value['weakly_ergodic_actually'])),
              output_filename='hist_A4')

    plot_hist(short_db_filename=f'{data_folder}simulations/smoothness/short_db_manycycles.pkl',
              filter_callable=lambda value: ((not value['single_use']) and (value['weakly_ergodic_actually'])),
              output_filename='hist_A5')

def make_hists_for_B_cases():
    plot_hist(short_db_filename=f'{data_folder}simulations/smoothness/short_db_B1-2.pkl',
              filter_callable=lambda value: (value['case'] == 'B1') and (not value['weakly_ergodic_actually']),
              output_filename='hist_B1')

    # THERE WERE NO INSTANCES WHERE THIS CASE OCCURS, SO IT IS OMITTED. CASE B1 WITH weakly_ergodic_actually IS NAMED
    # "B2" IN THE ARTICLE INSTEAD.
    # plot_hist(short_db_filename=f'{data_folder}simulations/smoothness/short_db_B1-2.pkl',
    #           filter_callable=lambda value: (value['case'] == 'B2') and (not value['weakly_ergodic_actually']),
    #           output_filename='hist_impossible')
    plot_hist(short_db_filename=f'{data_folder}simulations/smoothness/short_db_B1-2.pkl',
              filter_callable=lambda value: (value['case'] == 'B1') and (value['weakly_ergodic_actually']),
              output_filename='hist_B2')
    plot_hist(short_db_filename=f'{data_folder}simulations/smoothness/short_db_B1-5.pkl',
              filter_callable=lambda value: ((value['case'] == 'B3') and (value['weakly_ergodic_actually'])),
              output_filename='hist_B3')
    plot_hist(short_db_filename=f'{data_folder}simulations/smoothness/short_db_B1-5.pkl',
              filter_callable=lambda value: ((value['case'] == 'B4') and (value['weakly_ergodic_actually'])),
              output_filename='hist_B4')
    plot_hist(short_db_filename=f'{data_folder}simulations/smoothness/short_db_B1-5.pkl',
              filter_callable=lambda value: ((value['case'] == 'B5') and (value['weakly_ergodic_actually'])),
              output_filename='hist_B5')

if __name__ == '__main__':
    # analyze_one_reaction('reaction4ffc304e-0d60-4332-83d5-93f7f2bd0377')
    # make_hists_for_C_cases(case='C5')
    # make_hists_for_A_cases()
    make_hists_for_B_cases()