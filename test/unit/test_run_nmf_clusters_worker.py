import unittest
from unittest import TestCase
import numpy as np
import knpackage.toolbox as kn

import sample_clustering_toolbox_research_module as tstdata
import sample_clustering_toolbox as sctbx


class TestRun_nmf_clusters_worker(TestCase):
    def setUp(self):
        self.run_parameters = tstdata.get_test_paramters_dictionary()
        self.run_parameters["run_directory"] = '.'
        self.run_parameters["tmp_directory"] = kn.create_dir(self.run_parameters["run_directory"], 'tmp_cc_net_nmf')
        self.run_parameters["use_now_name"] = 0

    def tearDown(self):
        kn.remove_dir(self.run_parameters["tmp_directory"])
        del self.run_parameters

    def test_run_nmf_clusters_worker(self):
        self.run_parameters['number_of_bootstraps'] = 10
        self.run_parameters['method'] = self.run_parameters['method_1']

        cluster_rows = 3
        spreadsheet_consensus_mat = tstdata.get_square_3_cluster_spreadsheet(cluster_rows)
        spreadsheet_mat = tstdata.get_wide_3_cluster_spreadsheet(cluster_rows)

        for sample in range(0, int(self.run_parameters['number_of_bootstraps'])):
            sctbx.run_cc_nmf_clusters_worker(spreadsheet_mat, self.run_parameters, sample)

        linkage_matrix = np.zeros((spreadsheet_mat.shape[1], spreadsheet_mat.shape[1]))
        indicator_matrix = np.zeros((spreadsheet_mat.shape[1], spreadsheet_mat.shape[1]))
        consensus_matrix = sctbx.form_consensus_matrix(self.run_parameters, linkage_matrix, indicator_matrix)
        consensus_matrix[consensus_matrix != 0] = 1
        consensus_matrix = np.int_(consensus_matrix)

        cluster_difference = np.abs(spreadsheet_consensus_mat - consensus_matrix)
        cluster_difference = cluster_difference.sum()
        self.assertEqual(cluster_difference, 0, msg='nmf clustering failed')

if __name__ == '__main__':
    unittest.main()
