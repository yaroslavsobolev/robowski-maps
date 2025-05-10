import logging

from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import networkx as nx
import logging
logging.basicConfig(level=logging.ERROR)
import pickle
import uuid
import os
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map
# set highest protocol to 4
pickle.HIGHEST_PROTOCOL = 4

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

class ReactionNetwork:
    def __init__(self, N=10, single_use=True, weakly_ergodic=True, strict_backward=False, no_backward=False, min_k=0.001, max_k=1, random_seed=1337, no_cycles=False):
        self.g = nx.DiGraph()
        self.single_use = single_use
        self.weakly_ergodic = weakly_ergodic
        self.weakly_ergodic_actually = False
        self.strict_backward = strict_backward
        self.no_backward = no_backward
        self.no_cycles = no_cycles
        self.MIN_K = min_k
        self.MAX_K = max_k
        self.N = N
        self.g.add_nodes_from(range(N))
        self.molar_masses = np.ones(N)
        self.k_tensor = np.zeros((N, N, N))
        self.fist_order_K_matrix = np.zeros((N, N))
        self.reactions = list()
        self.set_of_used_reactants = set()
        self.condition_scans = []
        self.uuid_here = str(uuid.uuid4())
        self.random_seed = random_seed

    def add_reaction_of_second_order(self):
        if not self.weakly_ergodic:
            if not self.single_use:
                A_index = np.random.randint(0, self.N)
                B_index = np.random.randint(0, self.N)
            else:
                A_index = np.random.choice(list(set(range(self.N)) - self.set_of_used_reactants))
                B_index = np.random.choice(list(set(range(self.N)) - self.set_of_used_reactants))
        else:
            pairs_not_ergodically_connected = self.get_pairs_not_ergodically_connected()
            if self.single_use:
                pairs_not_ergodically_connected = [pair for pair in pairs_not_ergodically_connected if
                                                   ((pair[0] not in self.set_of_used_reactants) and (pair[
                                                       1] not in self.set_of_used_reactants))]
            if len(pairs_not_ergodically_connected) > 0:
                A_index, B_index = pairs_not_ergodically_connected[np.random.randint(0, len(pairs_not_ergodically_connected))]
            else:
                logging.info(f'All pairs are ergodically connected (single use : {self.single_use}). Choosing at random')
                if not self.single_use:
                    A_index = np.random.randint(0, self.N)
                    B_index = np.random.randint(0, self.N)
                else:
                    set_to_choose_from = list(set(range(self.N)) - self.set_of_used_reactants)
                    if len(set_to_choose_from) == 0:
                        logging.info(f'All reactants are used. Cannot add reaction.')
                        return True
                    else:
                        A_index = np.random.choice(list(set(range(self.N)) - self.set_of_used_reactants))
                        B_index = np.random.choice(list(set(range(self.N)) - self.set_of_used_reactants))

        set_of_possible_products = set(range(self.N)) - {A_index, B_index}
        all_ancestors = set(list(nx.ancestors(self.g, B_index)) + list(nx.ancestors(self.g, A_index)))
        if not self.no_backward:
            if not self.strict_backward:
                all_predecessors = set(list(self.g.predecessors(B_index)) + list(self.g.predecessors(A_index)))
            else:
                all_predecessors = set(list(self.g.predecessors(B_index))).intersection(set(list(self.g.predecessors(A_index))))
            ancestors_but_not_predecessors = all_ancestors - all_predecessors
            set_of_possible_products = set_of_possible_products - ancestors_but_not_predecessors
        else:
            set_of_possible_products = set_of_possible_products - all_ancestors

        if len(set_of_possible_products) == 0:
            logging.info(f'No possible products for {A_index} and {B_index}')
            return

        n_products = np.random.randint(1, len(set_of_possible_products) + 1)
        products = np.random.choice(list(set_of_possible_products), n_products, replace=False)
        local_ks = []
        for product in products:
            if not self.g.has_edge(A_index, product):
                self.g.add_edge(A_index, product)
            if not self.g.has_edge(B_index, product):
                self.g.add_edge(B_index, product)
            k_here = np.random.uniform(self.MIN_K, self.MAX_K)
            min_index = min(A_index, B_index)
            max_index = max(A_index, B_index)
            self.k_tensor[A_index, min_index, max_index] += -1 * k_here * self.molar_masses[A_index]
            self.k_tensor[B_index, min_index, max_index] += -1 * k_here * self.molar_masses[B_index]
            self.k_tensor[product, min_index, max_index] += k_here * (self.molar_masses[A_index] + self.molar_masses[B_index])
            logging.info(f'Added reaction {A_index} + {B_index} -> {product} with local k {k_here:.5f}')
            self.reactions.append({'substrates': [A_index, B_index], 'product': product, 'k': k_here})
            local_ks.append(k_here)

        self.set_of_used_reactants.add(A_index)
        self.set_of_used_reactants.add(B_index)


    def add_one_random_reversible_mixed_order_reaction(self):
        if not self.weakly_ergodic:
            if not self.single_use:
                A_index = np.random.randint(0, self.N)
                B_index = np.random.randint(0, self.N)
            else:
                set_to_choose_from = list(set(range(self.N)) - self.set_of_used_reactants)
                if len(set_to_choose_from) == 0:
                    logging.info(f'All reactants are used. Cannot add reaction.')
                    return True
                else:
                    A_index = np.random.choice(set_to_choose_from)
                    B_index = np.random.choice(set_to_choose_from)
        else:
            pairs_not_ergodically_connected = self.get_pairs_not_ergodically_connected()
            if self.single_use:
                pairs_not_ergodically_connected = [pair for pair in pairs_not_ergodically_connected if
                                                   ((pair[0] not in self.set_of_used_reactants) and (pair[
                                                       1] not in self.set_of_used_reactants))]
            if len(pairs_not_ergodically_connected) > 0:
                A_index, B_index = pairs_not_ergodically_connected[np.random.randint(0, len(pairs_not_ergodically_connected))]
            else:
                logging.info(f'All pairs are ergodically connected (single use : {self.single_use}). Choosing at random')
                if not self.single_use:
                    A_index = np.random.randint(0, self.N)
                    B_index = np.random.randint(0, self.N)
                else:
                    set_to_choose_from = list(set(range(self.N)) - self.set_of_used_reactants)
                    if len(set_to_choose_from) == 0:
                        logging.info(f'All reactants are used. Cannot add reaction.')
                        return True
                    else:
                        A_index = np.random.choice(set_to_choose_from)
                        B_index = np.random.choice(set_to_choose_from)

        logging.info(f'Adding reaction between {A_index} and {B_index}')
        set_of_possible_products = set(range(self.N)) - {A_index, B_index}
        logging.info(f'Ancestors of {A_index}: {list(nx.ancestors(self.g, A_index))}')
        logging.info(f'Ancestors of {B_index}: {list(nx.ancestors(self.g, B_index))}')
        all_ancestors = set(list(nx.ancestors(self.g, B_index)) + list(nx.ancestors(self.g, A_index)))
        logging.info(f'All ancestors: {all_ancestors}')
        if self.no_cycles:
            set_of_possible_products = set_of_possible_products - all_ancestors
        logging.info(f'Possible products: {set_of_possible_products}')

        if len(set_of_possible_products) == 0:
            logging.info(f'No possible products for {A_index} and {B_index}. Cannot add reaction.')
            return

        self.set_of_used_reactants.add(A_index)
        self.set_of_used_reactants.add(B_index)

        if len(set_of_possible_products) == 1:
            n_products = 1 # only one possible product left
        else:
            n_products = np.random.randint(1, 3) # random number of products, either 1 or 2

        if n_products == 2:
            # both forward and backward reactions are second-order
            products = np.random.choice(list(set_of_possible_products), n_products, replace=False)
            for product in products:
                if not self.g.has_edge(A_index, product):
                    self.g.add_edge(A_index, product)
                if not self.g.has_edge(B_index, product):
                    self.g.add_edge(B_index, product)
            k_here = np.random.uniform(self.MIN_K, self.MAX_K)
            min_index = min(A_index, B_index)
            max_index = max(A_index, B_index)
            self.k_tensor[A_index, min_index, max_index] += -1 * k_here * self.molar_masses[A_index]
            self.k_tensor[B_index, min_index, max_index] += -1 * k_here * self.molar_masses[B_index]
            # mass is spread evenly between the two products, if molar masses are equal to 1
            self.k_tensor[products[0], min_index, max_index] += 1 / 2 * k_here * (
                    self.molar_masses[A_index] + self.molar_masses[B_index])
            self.k_tensor[products[1], min_index, max_index] += 1 / 2 * k_here * (
                    self.molar_masses[A_index] + self.molar_masses[B_index])
            logging.info(f'Added reaction {A_index} + {B_index} -> {products} with local k {k_here:.5f}')
            self.reactions.append({'substrates': [A_index, B_index], 'product': products[0], 'k': k_here})
            self.reactions.append({'substrates': [A_index, B_index], 'product': products[1], 'k': k_here})

            # add a reverse second-order reaction between product 1 and product 2
            product_1, product_2 = products
            k_here = np.random.uniform(self.MIN_K, self.MAX_K)
            for product in products:
                if not self.g.has_edge(product, A_index):
                    self.g.add_edge(product, A_index)
                if not self.g.has_edge(product, B_index):
                    self.g.add_edge(product, B_index)
            min_index = min(product_1, product_2)
            max_index = max(product_1, product_2)
            self.k_tensor[product_1, min_index, max_index] += -1 * k_here * self.molar_masses[product_1]
            self.k_tensor[product_2, min_index, max_index] += -1 * k_here * self.molar_masses[product_2]
            self.k_tensor[A_index, min_index, max_index] += 1 / 2 * k_here * (
                    self.molar_masses[product_1] + self.molar_masses[product_2])
            self.k_tensor[B_index, min_index, max_index] += 1 / 2 * k_here * (
                    self.molar_masses[product_1] + self.molar_masses[product_2])
            self.reactions.append({'substrates': [product_1, product_2], 'product': A_index, 'k': k_here})
            self.reactions.append({'substrates': [product_1, product_2], 'product': B_index, 'k': k_here})
            self.set_of_used_reactants.add(products[0])
            self.set_of_used_reactants.add(products[1])

        if n_products == 1:
            # forward reaction is second-order, backward reaction is first-order
            product = np.random.choice(list(set_of_possible_products))
            if not self.g.has_edge(A_index, product):
                self.g.add_edge(A_index, product)
            if not self.g.has_edge(B_index, product):
                self.g.add_edge(B_index, product)
            k_here = np.random.uniform(self.MIN_K, self.MAX_K)
            min_index = min(A_index, B_index)
            max_index = max(A_index, B_index)
            self.k_tensor[A_index, min_index, max_index] += -1 * k_here * self.molar_masses[A_index]
            self.k_tensor[B_index, min_index, max_index] += -1 * k_here * self.molar_masses[B_index]
            self.k_tensor[product, min_index, max_index] += k_here * (self.molar_masses[A_index] + self.molar_masses[B_index])
            self.reactions.append({'substrates': [A_index, B_index], 'product': product, 'k': k_here
                                      })
            logging.info(f'Added reaction {A_index} + {B_index} -> {product} with local k {k_here:.5f}')

            # reverse reaction is first-order
            k_here = np.random.uniform(self.MIN_K, self.MAX_K)
            substrate_index = product
            if not self.g.has_edge(substrate_index, A_index):
                self.g.add_edge(substrate_index, A_index)
            if not self.g.has_edge(substrate_index, B_index):
                self.g.add_edge(substrate_index, B_index)
            self.fist_order_K_matrix[substrate_index, substrate_index] += -1 * k_here * self.molar_masses[
                substrate_index]
            self.fist_order_K_matrix[A_index, substrate_index] += 1/2 * k_here * self.molar_masses[substrate_index]
            self.fist_order_K_matrix[B_index, substrate_index] += 1/2 * k_here * self.molar_masses[substrate_index]
            logging.info(f'Added 1st order reaction {substrate_index} -> {A_index} + {B_index} with local k {k_here:.5f}')
            self.reactions.append({'substrates': [substrate_index], 'product': A_index, 'k': k_here})
            self.reactions.append({'substrates': [substrate_index], 'product': A_index, 'k': k_here})
            self.set_of_used_reactants.add(product)

        logging.info(f'Graph has number of edges: {self.g.number_of_edges()}')
        logging.info(f'List of all edges of the graph: {self.g.edges}')


    def add_specific_second_order_reaction(self, A_index, B_index, products):

        k_here = np.random.uniform(self.MIN_K, self.MAX_K)
        sum_of_product_molar_masses = 0
        for product in products:
            sum_of_product_molar_masses += self.molar_masses[product]

        min_index = min(A_index, B_index)
        max_index = max(A_index, B_index)
        self.k_tensor[A_index, min_index, max_index] += -1 * k_here * self.molar_masses[A_index]
        self.k_tensor[B_index, min_index, max_index] += -1 * k_here * self.molar_masses[B_index]

        for product in products:
            if not self.g.has_edge(A_index, product):
                self.g.add_edge(A_index, product)
            if not self.g.has_edge(B_index, product):
                self.g.add_edge(B_index, product)
            self.k_tensor[product, min_index, max_index] += k_here * (self.molar_masses[A_index] + self.molar_masses[B_index]) * self.molar_masses[product] / sum_of_product_molar_masses
            logging.info(f'Added reaction {A_index} + {B_index} -> {product} with local k {k_here:.5f}')
            self.reactions.append({'substrates': [A_index, B_index], 'product': product, 'k': k_here})

        self.set_of_used_reactants.add(A_index)
        self.set_of_used_reactants.add(B_index)


    def add_specific_first_order_reaction(self, substrate_index, product_index):

        k_here = np.random.uniform(self.MIN_K, self.MAX_K)
        if not self.g.has_edge(substrate_index, product_index):
            self.g.add_edge(substrate_index, product_index)
        self.fist_order_K_matrix[substrate_index, substrate_index] += -1 * k_here * self.molar_masses[substrate_index]
        self.fist_order_K_matrix[product_index, substrate_index] += 1 * k_here * self.molar_masses[substrate_index]
        logging.info(f'Added 1st order reaction {substrate_index} -> {product_index} with local k {k_here:.5f}')
        self.reactions.append({'substrates': [substrate_index], 'product': product_index, 'k': k_here})

        self.set_of_used_reactants.add(substrate_index)

    def get_pairs_not_ergodically_connected(self):
        pair_matrix = np.zeros((self.N, self.N))
        for node in self.g.nodes:
            ancestors_plus_node = list(nx.ancestors(self.g, node)) + [node]
            for ancestor1 in ancestors_plus_node:
                for ancestor2 in ancestors_plus_node:
                    pair_matrix[ancestor1, ancestor2] = 1
        return np.argwhere(pair_matrix == 0)

    def right_hand_side_of_ode_slow(self, t, y):
        # more transparent but very slow implementation of right_hand_side_of_ode(self, t, y)
        deriv_scale = 1
        derivs = np.zeros_like(y)
        for i in range(len(y)):
            for k in range(0, len(y)):
                for j in range(0, k + 1):
                    derivs[i] += self.k_tensor[i, j, k] * y[j] * y[k]
        return derivs * deriv_scale

    def right_hand_side_of_ode(self, t, y):
        return np.einsum('ijk,j,k->i', self.k_tensor, y, y) + (self.fist_order_K_matrix @ y)

    def solve_ode(self, t, y0, first_step=0.000001, max_step=100, atol=1e-8, rtol=1e-7):
        return solve_ivp(self.right_hand_side_of_ode, (0, t[-1]), y0, t_eval=t,
                         first_step=first_step, max_step=max_step, atol=atol, rtol=rtol)

    def add_reactions_of_second_order(self, n_reactions, ignore_cycles=False):
        for i in range(n_reactions):
            self.add_reaction_of_second_order()
        # normalize the tensor by the largest absolute value of components
        max_abs = np.max(np.abs(self.k_tensor))
        # print(f'Maximum absolute value of k tensor: {max_abs}')
        self.k_tensor = self.k_tensor / max_abs

        not_erodically_connected_pairs = self.get_pairs_not_ergodically_connected()
        if len(not_erodically_connected_pairs) > 0:
            self.weakly_ergodic_actually = False
            logging.info(f'Graph is not weakly ergodic')
        else:
            self.weakly_ergodic_actually = True
            logging.info(f'Graph is weakly ergodic')

        if not ignore_cycles:
            # maximum length of simple cycles
            simple_cycles = sorted(nx.simple_cycles(self.g))
            if len(simple_cycles) == 0:
                max_simple_cycle_length = 0
            else:
                max_simple_cycle_length = max([len(cycle) for cycle in simple_cycles])
            assert max_simple_cycle_length <= 2
            self.max_simple_cycle_length = max_simple_cycle_length


    def add_reactions_of_mixed_order(self, n_reactions, ignore_cycles=False):
        for i in range(n_reactions):
            self.add_one_random_reversible_mixed_order_reaction()
        # normalize the tensor by the largest absolute value of components
        max_abs = np.max(np.abs(self.k_tensor))
        # print(f'Maximum absolute value of k tensor: {max_abs}')
        maxk = max(np.max(np.abs(self.k_tensor)), np.max(np.abs(self.fist_order_K_matrix)))
        self.k_tensor = self.k_tensor / maxk
        self.fist_order_K_matrix = self.fist_order_K_matrix / maxk

        not_erodically_connected_pairs = self.get_pairs_not_ergodically_connected()
        if len(not_erodically_connected_pairs) > 0:
            self.weakly_ergodic_actually = False
            logging.info(f'Graph is not weakly ergodic')
        else:
            self.weakly_ergodic_actually = True
            logging.info(f'Graph is weakly ergodic')

        # if not ignore_cycles:
        #     # maximum length of simple cycles
        #     simple_cycles = sorted(nx.simple_cycles(self.g))
        #     if len(simple_cycles) == 0:
        #         max_simple_cycle_length = 0
        #     else:
        #         max_simple_cycle_length = max([len(cycle) for cycle in simple_cycles])
        #
        #     logging.info(f'Maximum simple cycle length: {max_simple_cycle_length}')
        #     if max_simple_cycle_length > 4:
        #         logging.info(f'Simple cycles: {simple_cycles}')
        #
        #     assert max_simple_cycle_length <= 4
        #
        #     self.max_simple_cycle_length = max_simple_cycle_length


    def scan_over_init_conditions(self, n_initial_conditions, min_c=0.01, delta_relative=0.001):
        # sample time logarithmically between 0.0001 and 1000
        ts = np.logspace(-3, 3.5, 200)

        # y0 = np.ones(N)
        max_over_t = []
        for initial_condition_index in range(n_initial_conditions):
            y0 = np.random.uniform(min_c, 1, self.N)
            y0 = y0 / np.sum(y0)
            solution = self.solve_ode(ts, y0)
            concentrations = solution.y.T
            times = solution.t

            # sum of concentration at the end of the simulation
            # sum_concentrations = np.sum(concentrations[-1])
            # print(f'Sum of concentrations at the end of the simulation: {sum_concentrations}')

            # plt.plot(solution.t, solution.y.T, label=range(N))
            # # make log scale in x axis
            # plt.xscale('log')
            # plt.xlabel('Time')
            # plt.ylabel('Concentration')
            # plt.legend()
            # plt.show()

            # evaluating partial derivatives of elements of y at the end of the simulation by the initial concentrations y0
            partial_derivs = []
            for i in range(self.N):
                y0_plus_delta = y0.copy()
                delta = y0_plus_delta[i] * delta_relative
                y0_plus_delta[i] += delta
                solution_plus_delta = self.solve_ode(ts, y0_plus_delta)
                y_deriv = (solution_plus_delta.y.T - solution.y.T) / delta
                partial_derivs.append(y_deriv.copy())
                # plt.plot(solution.t, y_deriv, label=i)

            partial_derivs = np.array(partial_derivs)
            # highest partial derivative for each given time
            max_partial_derivs = np.max(np.abs(partial_derivs), axis=(0, 2))
            max_over_t.append(max_partial_derivs)
            # self.condition_scans.append({'y0': y0, 'partial_derivs': partial_derivs, 'max_partial_derivs_over_t': max_partial_derivs})

        # make a max over all initial conditions for max_over_t
        max_over_t = np.array(max_over_t)
        max_over_t = np.max(max_over_t, axis=0)

        self.max_over_t = np.copy(max_over_t)
        # print(f'Maximum partial derivative: {np.max(np.abs(partial_derivs))}')
        # print(f'It is achieved for index {np.unravel_index(np.argmax(np.abs(partial_derivs)), partial_derivs.shape)}')

    def scan_over_special_init_conditions(self, n_initial_conditions, substrate_ids, save_condition_with_highest_deriv=False, verbose=False, custom_y0_list=None,
                                          min_c=0.01):
        # sample time logarithmically between 0.0001 and 1000
        ts = np.logspace(-3, 3.5, 200)

        # y0 = np.ones(N)
        global_max_deriv = 0
        max_over_t = []
        for initial_condition_index in range(n_initial_conditions):
            if verbose:
                print(f'Initial condition {initial_condition_index} of {n_initial_conditions}')
            if custom_y0_list is not None:
                y0 = custom_y0_list[initial_condition_index]
            else:
                y0 = np.random.uniform(min_c, 1, self.N)
                # set all values that are not in substrate_ids list to zero
                for i in range(self.N):
                    if i not in substrate_ids:
                        y0[i] = 0
                y0 = y0 / np.sum(y0)

            solution = self.solve_ode(ts, y0)
            concentrations = solution.y.T
            times = solution.t
            # evaluating partial derivatives of elements of y at the end of the simulation by the initial concentrations y0
            partial_derivs = []
            delta_relative = 0.001
            for i in substrate_ids:
                y0_plus_delta = y0.copy()
                delta = y0_plus_delta[i] * delta_relative
                y0_plus_delta[i] += delta
                solution_plus_delta = self.solve_ode(ts, y0_plus_delta)
                y_deriv = (solution_plus_delta.y.T - solution.y.T) / delta
                partial_derivs.append(y_deriv.copy())
                # plt.plot(solution.t, y_deriv, label=i)

            partial_derivs = np.array(partial_derivs)
            # highest partial derivative for each given time
            max_partial_derivs = np.max(np.abs(partial_derivs), axis=(0, 2))
            max_over_t.append(max_partial_derivs)

            global_max_here = np.max(np.max(max_partial_derivs))
            if global_max_here > global_max_deriv:
                global_max_deriv = global_max_here
                if save_condition_with_highest_deriv:
                    if verbose:
                        print(f'Found new global maximum: {global_max_here}')
                    self.condition_with_highest_deriv = {'y0': y0, 'max_partial_derivs_over_t': max_partial_derivs}
            # self.condition_scans.append({'y0': y0, 'partial_derivs': partial_derivs, 'max_partial_derivs_over_t': max_partial_derivs})

        # make a max over all initial conditions for max_over_t
        max_over_t = np.max(np.array(max_over_t), axis=0)

        self.max_over_t = np.copy(max_over_t)
        self.max_overall = np.max(np.max(max_over_t))

    def pickle_results(self, target_path=None):
        if target_path is None:
            folder_name = f'simulations/reaction{self.uuid_here}'
        else:
            folder_name = target_path
        # make folder if not exists
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)
        # pickle the rs to this folder
        with open(folder_name + '/reaction_network.pkl', 'wb') as f:
            pickle.dump(self, f)


def construct_short_db():
    short_db = {}
    # target_folder = f'{data_folder}simulations/smoothness/'
    target_folder = 'simulations/'
    # find all folders starting with 'reaction' in the target folder uding glob
    import glob
    folders = glob.glob(target_folder + 'reaction*')
    print(f'Found {len(folders)} folders')
    # load all reaction networks from these folders
    max_d = 0
    for folder in tqdm(folders):
        folder_name_without_path = folder.split('/')[-1]
        rs2 = pickle.load(open(folder + '/reaction_network.pkl', 'rb'))
        # maximum length of simple cycles
        # simple_cycles = sorted(nx.simple_cycles(rs2.g))
        # if len(simple_cycles) == 0:
        #     max_simple_cycle_length = 0
        # else:
        #     max_simple_cycle_length = max([len(cycle) for cycle in simple_cycles])
        # print(f'Found {len(simple_cycles)} simple cycles, the longest one has length {max_simple_cycle_length}')
        max_d = max(max_d, np.max(rs2.max_over_t))
        # print(f'Maxoverall: {rs2.max_overall}')
        # print(f'max_over_t: {rs2.max_over_t}')

        short_db[folder_name_without_path] = {'max_over_t': rs2.max_over_t, 'reactions': rs2.reactions,
                                              'k_tensor': rs2.k_tensor, 'single_use': rs2.single_use,
                                              'weakly_ergodic': rs2.weakly_ergodic,
                                              'strict_backward': rs2.strict_backward,
                                              'weakly_ergodic_actually': rs2.weakly_ergodic_actually,
                                              'no_backward': rs2.no_backward}

        # if rs2.case object exists
        if hasattr(rs2, 'case'):
            short_db[folder_name_without_path]['case'] = rs2.case
        if hasattr(rs2, 'max_simple_cycle_length'):
            short_db[folder_name_without_path]['max_simple_cycle_length'] = rs2.max_simple_cycle_length
        if hasattr(rs2, 'no_cycles'):
            short_db[folder_name_without_path]['no_cycles'] = rs2.no_cycles

        # 'max_simple_cycle_length': max_simple_cycle_length}

    print(f'Maximum value of max_over_t: {max_d}')

    # pickle the short_db
    with open('short_db.pkl', 'wb') as f:
        pickle.dump(short_db, f)


def load_short_db(suffix=''):
    with open(f'{data_folder}simulations/smoothness/short_db{suffix}.pkl', 'rb') as f:
        short_db = pickle.load(f)
    return short_db


if __name__ == '__main__':
    target_folder = f'{data_folder}simulations/smoothness/old_simulations_with_C_from_0_to_1'
    with open(f'{data_folder}simulations/smoothness/short_db_noback.pkl', 'rb') as f:
        short_db = pickle.load(f)

    # plot the max_over_t for all reactions that are weakly ergodic, strict backward and single use
    ts = np.logspace(-3, 3.5, 200)
    max_over_ts = []
    for key, value in short_db.items():
        if (value['single_use']):
            max_over_ts.append(value['max_over_t'])
            # if np.max(value['max_over_t']) > 7.83:
            #     print(key)
    print(f'Found {len(max_over_ts)} reactions that match the criteria')
    max_over_ts = np.array(max_over_ts)
    max_over_ts_global = np.max(max_over_ts, axis=0)
    max_for_each_reaction = np.max(max_over_ts, axis=1)
    f1 = plt.figure(figsize=(5, 1.5))
    plt.hist(max_for_each_reaction, bins=100, density=True, color='grey')
    plt.axvline(np.max(max_over_ts_global), color='red')
    plt.xlabel('Highest partial derivative for a given reaction')
    plt.ylabel('Probability density')
    plt.tight_layout()
    plt.show()
    print(f'Maximum value of max_over_ts_global: {np.max(max_over_ts_global)}')
    np.save(f'{target_folder}/max_over_ts_global.npy', max_over_ts_global)
    plt.scatter(ts, max_over_ts_global)
    plt.xscale('log')
    # plt.axhline(1, color='red')
    # plt.yscale('log')
    plt.xlabel('Reaction time, a.u. (note the logarithmic scale)')
    plt.ylabel('Maximum absolute partial derivative\n(maximum across all reactions)')
    # plot the exp(4*ts) curve
    # plt.plot(ts, np.exp(4 * ts), color='black')
    plt.show()

    # # # make a 2d heatmap histogram for the max_over_ts, with axes corresponding to max_over_ts axes
    # # for i in range(max_over_ts.shape[0]):
    # #     plt.scatter(ts, max_over_ts[i, :], alpha=0.1)
    # # plt.xscale('log')
    # # plt.show()

    interesting_reaction = 'reactioncd4d40a4-ceae-4fee-bef6-9a26b06c12bc'
    threshold = 7.83
    y0 = np.zeros(10)
    re2 = pickle.load(open(f'{target_folder}/{interesting_reaction}/reaction_network.pkl', 'rb'))
    for condition_id, condition_scan in enumerate(re2.condition_scans):
        if np.max(condition_scan['max_partial_derivs_over_t']) > threshold:
            print(f'Found condition with max partial derivative over t > {threshold}, the index is {condition_id}')
            y0 = condition_scan['y0']
    print(y0)

    # Plot the graph g
    # maximum length of simple cycles
    simple_cycles = sorted(nx.simple_cycles(re2.g))
    if len(simple_cycles) == 0:
        max_simple_cycle_length = 0
    else:
        max_simple_cycle_length = max([len(cycle) for cycle in simple_cycles])
    print(f'Found {len(simple_cycles)} simple cycles, the longest one has length {max_simple_cycle_length}')

    # list weakly connected components
    weakly_connected_components = list(nx.weakly_connected_components(re2.g))
    print(f'Found {len(weakly_connected_components)} weakly connected components')
    for i, component in enumerate(weakly_connected_components):
        print(f'Component {i} is {component}')

    # strongly connected components:
    strongly_connected_components = list(nx.strongly_connected_components(re2.g))
    print(f'Found {len(strongly_connected_components)} strongly connected components')
    for i, component in enumerate(strongly_connected_components):
        print(f'Component {i} is {component}')

    # print(re2.g.edges())
    nx.draw(re2.g, with_labels=True, pos=nx.spring_layout(re2.g))
    plt.show()

    # remove node 6 from the graph
    re2.g.remove_node(6)
    print('Pairs not ergodically connected, besides node 6:')
    print(re2.get_pairs_not_ergodically_connected())
    # print(re2.g.edges())
    nx.draw(re2.g, with_labels=True, pos=nx.spring_layout(re2.g))
    plt.show()

    for node in re2.g.nodes:
        re2.g.nodes[node]['URL'] = str(node)

    nx.write_graphml(re2.g, 'graph.graphml')

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

    # re2 = ReactionNetwork(N=10, single_use=True, weakly_ergodic=True, strict_backward=True, min_k=0.001, max_k=1)
    # re2.k_tensor = short_db[interesting_reaction]['k_tensor']
    # re2.