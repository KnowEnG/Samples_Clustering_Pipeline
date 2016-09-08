import unittest
from unittest import TestCase
import numpy as np
import knpackage.toolbox as kn

import sample_clustering_toolbox_research_module as tstdata
import sample_clustering_toolbox as sctbx


class TestRun_nmf_clusters_worker(TestCase):

    def test_run_nmf_clusters_worker(self):
        tmp_dir = 'tmp_cc_net_nmf'
        run_parameters = tstdata.get_test_paramters_dictionary()
        run_parameters['number_of_bootstraps'] = 10
        run_parameters["run_directory"] = '.'
        run_parameters["tmp_directory"] = kn.create_dir(run_parameters["run_directory"], tmp_dir)
        run_parameters['method'] = run_parameters['method_1']
        run_parameters['display_clusters'] = 0
        run_parameters["use_now_name"] = 0

        cluster_rows = 5
        spreadsheet_consensus_mat = tstdata.get_square_3_cluster_spreadsheet(cluster_rows)
        spreadsheet_mat = tstdata.get_wide_3_cluster_spreadsheet(cluster_rows)

        for sample in range(0, int(run_parameters['number_of_bootstraps'])):
            sctbx.run_nmf_clusters_worker(spreadsheet_mat, run_parameters, sample)

        linkage_matrix = np.zeros((spreadsheet_mat.shape[1], spreadsheet_mat.shape[1]))
        indicator_matrix = np.zeros((spreadsheet_mat.shape[1], spreadsheet_mat.shape[1]))
        consensus_matrix = sctbx.form_consensus_matrix(run_parameters, linkage_matrix, indicator_matrix)
        consensus_matrix[consensus_matrix != 0] = 1
        consensus_matrix = np.int_(consensus_matrix)

        kn.remove_dir(run_parameters["tmp_directory"])

        cluster_difference = np.abs(spreadsheet_consensus_mat - consensus_matrix)
        cluster_difference = cluster_difference.sum()
        self.assertEqual(cluster_difference, 0, msg='nmf clustering failed')

if __name__ == '__main__':
    unittest.main()

"""
def run_nmf_clusters_worker(spreadsheet_mat, run_parameters, sample):
    #""Worker to execute nmf_clusters in a single process

    Args:
        spreadsheet_mat: genes x samples matrix.
        run_parameters: dictionary of run-time parameters.
        sample: each loops.

    Returns:
        None

    #""
    sample_random, sample_permutation = kn.sample_a_matrix(
        spreadsheet_mat, np.float64(run_parameters["rows_sampling_fraction"]),
        np.float64(run_parameters["cols_sampling_fraction"]))

    h_mat = kn.perform_nmf(sample_random, run_parameters)
    save_a_clustering_to_tmp(h_mat, sample_permutation, run_parameters, sample)

    print('bootstrap {} of {}'.format(sample + 1, run_parameters["number_of_bootstraps"]))
"""
