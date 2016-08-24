# -*- coding: utf-8 -*-
"""
Created on WEd Aug  17 2016

@author: del
@author: The Gene Sets Characterization dev team
"""

import unittest
import numpy as np
#import scipy.sparse as spar
#import scipy.stats as stats

import sys
sample_clustering_unit_test_directory = '/Users/lanier4/PycharmProjects/Samples_Clustering_Pipeline/src'
sys.path.extend(sample_clustering_unit_test_directory)
import sample_clustering_toolbox as sctbx

def synthesize_random_network(network_dim, n_nodes):
    """ symmetric random adjacency matrix from random set of nodes
    Args:
        network_dim: number of rows and columns in the symmetric output matrix
        n_nodes: number of connections (approximate because duplicates are ignored)
    Returns:
        network: a symmetric adjacency matrix (0 or 1 in network_dim x network_dim matrix)
    """
    network = np.zeros((network_dim, network_dim))
    col_0 = np.random.randint(0, network_dim, n_nodes)
    col_1 = np.random.randint(0, network_dim, n_nodes)
    for node in range(0, n_nodes):
        if col_0[node] != col_1[node]:
            network[col_0[node], col_1[node]] = 1
    network = network + network.T
    network[network != 0] = 1

    return network

def get_cluster_indices_list(a_arr):
    """ get the list of sets of positive integers in the input array where a set
        is the index of where equal values occur for all equal values in the array

    Args:
        a_arr: array of positive integers

    Returns:
        cluster_list: list of lists where each list is the indecies of the members
            of that set, and the lists are ordered by the first member of each.
    """
    idx_arr = np.arange(0, a_arr.size)
    a_arr_unique = np.unique(a_arr)
    tmp_list = []
    for v in a_arr_unique:
        tmp_list.append(idx_arr[a_arr == v])

    len_tmp_list = len(tmp_list)
    first_member_array = np.int_(np.zeros(len_tmp_list))
    for m in range(0, len_tmp_list):
        tmp = tmp_list[m]
        first_member_array[m] = int(tmp[0])

    list_order = np.int_(np.argsort(first_member_array))
    cluster_list = []
    for t in list_order:
        cluster_list.append(tmp_list[t])

    return cluster_list

def sets_a_eq_b(a, b):
    """ check that all indices of equal values in a
        are same sets as indices of equal values in b
    Args:
        a: array of cluster assignments
        b: array of cluster assignments - same size or will return false
    Returns:
        True or False: array a indices of equal value
            are the same as array b indices of equal values
    """
    a_u = np.unique(a)
    b_u = np.unique(b)
    if len(a_u) != len(b_u):
        return False
    else:
        a_list = get_cluster_indices_list(a)
        b_list = get_cluster_indices_list(b)
        if len(b) != len(a):
            return False
        else:
            n_here = 0
            for a_set in a_list:
                if (len(a_set) != len(b_list[n_here])):
                    return False
                elif sum(np.int_(a_set != b_list[n_here])) != 0:
                    return False
                else:
                    n_here += 1
    return True


class sample_clustering_toolbox_test(unittest.TestCase):
    def get_run_parameters(self):
        run_parameters = {'test_directory': '/Users/lanier4/BigDataTank/nbs_run',
                          'k': 3, 'number_of_iteriations_in_rwr': 100,
                          'obj_fcn_chk_freq': 50,
                          'it_max': 10000,
                          'h_clust_eq_limit': 100,
                          'restart_tolerance': 0.0001,
                          'lmbda': 1400,
                          'percent_sample': 0.8,
                          'number_of_bootstraps': 3,
                          'display_clusters': 1,
                          'restart_probability': 0.7,
                          'verbose': 1,
                          'use_now_name': 1000000}

        return run_parameters

    def test_timestamp_filename(self):
        """ test sample_clustering_toolbox function timestamp_filename input switch
        """
        run_parameters = self.get_run_parameters()
        run_parameters['use_now_name'] = 1
        name_base = 'base_name'
        name_extension = 'ext_name'
        f_name = sctbx.timestamp_filename(name_base, name_extension, run_parameters)
        f_nameII = sctbx.timestamp_filename(name_base, name_extension)
        self.assertNotEqual(f_name, f_nameII, msg='{} != {}'.format(f_name, f_nameII))


def suite():
    test_suite = unittest.TestSuite()
    test_suite.addTest(unittest.makeSuite(sample_clustering_toolbox_test))

    return test_suite
