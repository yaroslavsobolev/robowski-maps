from robowski.settings import *
from smoothness_v2 import *

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

def run_one_simulation(random_seed):
    np.random.seed(random_seed)
    rs = ReactionNetwork(N=10, single_use=True, weakly_ergodic=True, strict_backward=True, no_backward=True, min_k=0.001, max_k=1)
    rs.add_reactions_of_second_order(20)
    rs.scan_over_init_conditions(100)
    rs.pickle_results()

def run_one_simulation_2(random_seed):
    np.random.seed(random_seed)
    rs = ReactionNetwork(N=10, single_use=False, weakly_ergodic=True, strict_backward=True, no_backward=True, min_k=0.001, max_k=1)
    rs.add_reactions_of_second_order(20)
    rs.scan_over_init_conditions(100)
    rs.pickle_results()

def run_one_simulation_3(random_seed):
    np.random.seed(random_seed)
    rs = ReactionNetwork(N=10, single_use=True, weakly_ergodic=True, strict_backward=False, min_k=0.001, max_k=1)
    rs.add_reactions_of_second_order(20, ignore_cycles=True)
    rs.scan_over_init_conditions(100)
    rs.pickle_results()

def run_one_simulation_4(random_seed):
    np.random.seed(random_seed)
    rs = ReactionNetwork(N=10, single_use=False, weakly_ergodic=True, strict_backward=False, min_k=0.001, max_k=1)
    rs.add_reactions_of_second_order(20, ignore_cycles=True)
    rs.scan_over_init_conditions(100)
    rs.pickle_results()


def run_one_simulation_B1(random_seed):
    np.random.seed(random_seed)
    rs = ReactionNetwork(N=10, single_use=True, weakly_ergodic=False, min_k=0.001, max_k=1, no_cycles=True)
    rs.add_reactions_of_mixed_order(20, ignore_cycles=True)
    rs.scan_over_init_conditions(100)
    rs.case = 'B1'
    rs.pickle_results()

def run_one_simulation_B2(random_seed):
    np.random.seed(random_seed)
    rs = ReactionNetwork(N=10, single_use=False, weakly_ergodic=False, min_k=0.001, max_k=1, no_cycles=True)
    rs.add_reactions_of_mixed_order(20, ignore_cycles=True)
    rs.scan_over_init_conditions(100)
    rs.case = 'B2'
    rs.pickle_results()

def run_one_simulation_B3(random_seed):
    np.random.seed(random_seed)
    rs = ReactionNetwork(N=10, single_use=False, weakly_ergodic=True, min_k=0.001, max_k=1, no_cycles=True)
    rs.add_reactions_of_mixed_order(20, ignore_cycles=True)
    rs.scan_over_init_conditions(100)
    rs.case = 'B3'
    rs.pickle_results()

def run_one_simulation_B4(random_seed):
    np.random.seed(random_seed)
    rs = ReactionNetwork(N=10, single_use=True, weakly_ergodic=True, min_k=0.001, max_k=1, no_cycles=False)
    rs.add_reactions_of_mixed_order(20, ignore_cycles=True)
    rs.scan_over_init_conditions(100)
    rs.case = 'B4'
    rs.pickle_results()

def run_one_simulation_B5(random_seed):
    np.random.seed(random_seed)
    rs = ReactionNetwork(N=10, single_use=False, weakly_ergodic=True, min_k=0.001, max_k=1, no_cycles=False)
    rs.add_reactions_of_mixed_order(20, ignore_cycles=True)
    rs.scan_over_init_conditions(100)
    rs.case = 'B5'
    rs.pickle_results()


def run_one_simulation_5(random_seed):
    rs = ReactionNetwork(N=14, single_use=False, weakly_ergodic=False, strict_backward=False, min_k=0.001, max_k=1, random_seed=random_seed)
    np.random.seed(random_seed)
    rs.add_specific_second_order_reaction(0, 1, [2, 3])
    rs.add_specific_second_order_reaction(2, 3, [0, 1])
    rs.add_specific_second_order_reaction(3, 4, [5, 6])
    rs.add_specific_second_order_reaction(5, 6, [3, 4])
    rs.add_specific_first_order_reaction(0, 7)
    rs.add_specific_first_order_reaction(7, 0)
    rs.add_specific_first_order_reaction(1, 8)
    rs.add_specific_first_order_reaction(8, 1)
    rs.add_specific_first_order_reaction(2, 9)
    rs.add_specific_first_order_reaction(9, 2)
    rs.add_specific_first_order_reaction(3, 10)
    rs.add_specific_first_order_reaction(10, 3)
    rs.add_specific_first_order_reaction(4, 11)
    rs.add_specific_first_order_reaction(11, 4)
    rs.add_specific_first_order_reaction(5, 12)
    rs.add_specific_first_order_reaction(12, 5)
    rs.add_specific_first_order_reaction(6, 13)
    rs.add_specific_first_order_reaction(13, 6)
    maxk = max(np.max(np.abs(rs.k_tensor)), np.max(np.abs(rs.fist_order_K_matrix)))
    rs.k_tensor = rs.k_tensor / maxk
    rs.fist_order_K_matrix = rs.fist_order_K_matrix / maxk
    rs.weakly_ergodic_actually = False
    rs.scan_over_special_init_conditions(100, [0, 1, 4], save_condition_with_highest_deriv=True)
    rs.case = 'C1'
    rs.max_simple_cycle_length = 0
    rs.pickle_results()


def run_one_simulation_6(random_seed):
    rs = ReactionNetwork(N=14, single_use=False, weakly_ergodic=False, strict_backward=False, min_k=0.001, max_k=1, random_seed=random_seed)
    np.random.seed(random_seed)
    rs.add_specific_second_order_reaction(0, 1, [2, 3])
    rs.add_specific_second_order_reaction(2, 3, [0, 1])
    rs.add_specific_second_order_reaction(3, 4, [5, 6])
    rs.add_specific_second_order_reaction(5, 6, [3, 4])
    rs.add_specific_first_order_reaction(0, 7)
    rs.add_specific_first_order_reaction(7, 0)
    rs.add_specific_first_order_reaction(1, 8)
    rs.add_specific_first_order_reaction(8, 1)
    rs.add_specific_first_order_reaction(2, 9)
    rs.add_specific_first_order_reaction(9, 2)
    rs.add_specific_first_order_reaction(3, 10)
    rs.add_specific_first_order_reaction(10, 3)
    rs.add_specific_first_order_reaction(4, 11)
    rs.add_specific_first_order_reaction(11, 4)
    rs.add_specific_first_order_reaction(5, 12)
    rs.add_specific_first_order_reaction(12, 5)
    rs.add_specific_first_order_reaction(6, 13)
    rs.add_specific_first_order_reaction(13, 6)
    rs.fist_order_K_matrix = rs.fist_order_K_matrix / 10
    maxk = max(np.max(np.abs(rs.k_tensor)), np.max(np.abs(rs.fist_order_K_matrix)))
    rs.k_tensor = rs.k_tensor / maxk
    rs.fist_order_K_matrix = rs.fist_order_K_matrix / maxk
    rs.weakly_ergodic_actually = False
    rs.scan_over_special_init_conditions(100, [0, 1, 4], save_condition_with_highest_deriv=True)
    rs.case = 'C2'
    rs.max_simple_cycle_length = 0
    rs.pickle_results()


def run_one_simulation_7(random_seed):
    rs = ReactionNetwork(N=14, single_use=False, weakly_ergodic=False, strict_backward=False, min_k=0.001, max_k=1, random_seed=random_seed)
    np.random.seed(random_seed)
    rs.add_specific_second_order_reaction(0, 1, [2, 3])
    rs.add_specific_second_order_reaction(2, 3, [0, 1])
    rs.add_specific_second_order_reaction(3, 4, [5, 6])
    rs.add_specific_second_order_reaction(5, 6, [3, 4])
    rs.add_specific_first_order_reaction(0, 7)
    rs.add_specific_first_order_reaction(7, 0)
    rs.add_specific_first_order_reaction(1, 8)
    rs.add_specific_first_order_reaction(8, 1)
    rs.add_specific_first_order_reaction(2, 9)
    rs.add_specific_first_order_reaction(9, 2)
    rs.add_specific_first_order_reaction(3, 10)
    rs.add_specific_first_order_reaction(10, 3)
    rs.add_specific_first_order_reaction(4, 11)
    rs.add_specific_first_order_reaction(11, 4)
    rs.add_specific_first_order_reaction(5, 12)
    rs.add_specific_first_order_reaction(12, 5)
    rs.add_specific_first_order_reaction(6, 13)
    rs.add_specific_first_order_reaction(13, 6)
    rs.fist_order_K_matrix = rs.fist_order_K_matrix / 100
    maxk = max(np.max(np.abs(rs.k_tensor)), np.max(np.abs(rs.fist_order_K_matrix)))
    rs.k_tensor = rs.k_tensor / maxk
    rs.fist_order_K_matrix = rs.fist_order_K_matrix / maxk
    rs.weakly_ergodic_actually = False
    rs.scan_over_special_init_conditions(100, [0, 1, 4], save_condition_with_highest_deriv=True)
    rs.case = 'C3'
    rs.max_simple_cycle_length = 0
    rs.pickle_results()

def run_one_simulation_8(random_seed):
    rs = ReactionNetwork(N=16, single_use=False, weakly_ergodic=False, strict_backward=False, min_k=0.001, max_k=1, random_seed=random_seed)
    np.random.seed(random_seed)
    rs.add_specific_second_order_reaction(0, 1, [2, 3])
    rs.add_specific_second_order_reaction(2, 3, [0, 1])

    rs.add_specific_first_order_reaction(3, 4)
    rs.add_specific_first_order_reaction(4, 3)

    rs.add_specific_second_order_reaction(4, 5, [6, 7])
    rs.add_specific_second_order_reaction(6, 7, [4, 5])

    # side reactions
    rs.add_specific_first_order_reaction(0, 8)
    rs.add_specific_first_order_reaction(8, 0)

    rs.add_specific_first_order_reaction(1, 9)
    rs.add_specific_first_order_reaction(9, 1)

    rs.add_specific_first_order_reaction(2, 10)
    rs.add_specific_first_order_reaction(10, 2)

    rs.add_specific_first_order_reaction(3, 11)
    rs.add_specific_first_order_reaction(11, 3)

    rs.add_specific_first_order_reaction(4, 12)
    rs.add_specific_first_order_reaction(12, 4)

    rs.add_specific_first_order_reaction(5, 13)
    rs.add_specific_first_order_reaction(13, 5)

    rs.add_specific_first_order_reaction(6, 14)
    rs.add_specific_first_order_reaction(14, 6)

    rs.add_specific_first_order_reaction(7, 15)
    rs.add_specific_first_order_reaction(15, 7)

    # rs.fist_order_K_matrix = rs.fist_order_K_matrix
    maxk = max(np.max(np.abs(rs.k_tensor)), np.max(np.abs(rs.fist_order_K_matrix)))
    rs.k_tensor = rs.k_tensor / maxk
    rs.fist_order_K_matrix = rs.fist_order_K_matrix / maxk
    rs.weakly_ergodic_actually = False
    rs.scan_over_special_init_conditions(100, [0, 1, 5], save_condition_with_highest_deriv=True)
    rs.case = 'C4'
    rs.max_simple_cycle_length = 0
    rs.pickle_results()


def run_one_simulation_9(random_seed):
    # CASE C5
    rs = ReactionNetwork(N=19, single_use=False, weakly_ergodic=False, strict_backward=False, min_k=0.001, max_k=1, random_seed=random_seed)
    np.random.seed(random_seed)
    rs.add_specific_second_order_reaction(0, 1, [2, 3])
    rs.add_specific_second_order_reaction(2, 3, [0, 1])
    rs.add_specific_second_order_reaction(3, 4, [5, 6])
    rs.add_specific_second_order_reaction(5, 6, [3, 4])
    rs.add_specific_first_order_reaction(0, 7)
    rs.add_specific_first_order_reaction(7, 0)
    rs.add_specific_first_order_reaction(1, 8)
    rs.add_specific_first_order_reaction(8, 1)
    rs.add_specific_first_order_reaction(2, 9)
    rs.add_specific_first_order_reaction(9, 2)
    rs.add_specific_first_order_reaction(3, 10)
    rs.add_specific_first_order_reaction(10, 3)
    rs.add_specific_first_order_reaction(4, 11)
    rs.add_specific_first_order_reaction(11, 4)
    rs.add_specific_first_order_reaction(5, 12)
    rs.add_specific_first_order_reaction(12, 5)
    rs.add_specific_first_order_reaction(6, 13)
    rs.add_specific_first_order_reaction(13, 6)

    rs.add_specific_first_order_reaction(13, 14)
    rs.add_specific_first_order_reaction(14, 13)

    rs.add_specific_first_order_reaction(14, 15)
    rs.add_specific_first_order_reaction(15, 14)

    rs.add_specific_first_order_reaction(13, 16)
    rs.add_specific_first_order_reaction(16, 13)

    rs.add_specific_first_order_reaction(14, 17)
    rs.add_specific_first_order_reaction(17, 14)

    rs.add_specific_first_order_reaction(15, 18)
    rs.add_specific_first_order_reaction(18, 15)

    maxk = max(np.max(np.abs(rs.k_tensor)), np.max(np.abs(rs.fist_order_K_matrix)))
    rs.k_tensor = rs.k_tensor / maxk
    rs.fist_order_K_matrix = rs.fist_order_K_matrix / maxk
    rs.weakly_ergodic_actually = False
    rs.scan_over_special_init_conditions(100, [0, 1, 4], save_condition_with_highest_deriv=True)
    rs.case = 'C5'
    rs.max_simple_cycle_length = 0
    rs.pickle_results()


def screen_one_reaction_over_conditions(interesting_reaction, target_folder = f'{data_folder}simulations/smoothness'):
    rs2 = pickle.load(open(f'{target_folder}/{interesting_reaction}/reaction_network.pkl', 'rb'))
    rs2.scan_over_special_init_conditions(1000, [0, 1, 4], save_condition_with_highest_deriv=True, verbose=True)
    print(rs2.condition_with_highest_deriv['y0'])
    print(rs2.max_overall)
    rs2.pickle_results(target_path=f'screened/C1_max_{interesting_reaction}')


if __name__ == '__main__':
    # r = process_map(run_one_simulation, range(10000), max_workers=77, chunksize=1)
    # r = process_map(run_one_simulation_2, range(10000), max_workers=77, chunksize=1)
    # r = process_map(run_one_simulation_3, range(10000), max_workers=77, chunksize=1)
    # r = process_map(run_one_simulation_4, range(10000), max_workers=77, chunksize=1)
    # r = process_map(run_one_simulation_6, range(10000), max_workers=77, chunksize=1)
    # r = process_map(run_one_simulation_7, range(10000), max_workers=77, chunksize=1)
    # r = process_map(run_one_simulation_9, range(10000), max_workers=77, chunksize=1)
    # r = process_map(run_one_simulation_9, range(10000), max_workers=77, chunksize=1)
    # screen_one_reaction_over_conditions(interesting_reaction='reaction4ffc304e-0d60-4332-83d5-93f7f2bd0377', target_folder='simulations_C1_C2')

    r = process_map(run_one_simulation_B1, range(10000), max_workers=77, chunksize=1)
    r = process_map(run_one_simulation_B2, range(10000), max_workers=77, chunksize=1)
    r = process_map(run_one_simulation_B3, range(10000), max_workers=77, chunksize=1)
    r = process_map(run_one_simulation_B4, range(10000), max_workers=77, chunksize=1)
    r = process_map(run_one_simulation_B5, range(10000), max_workers=77, chunksize=1)