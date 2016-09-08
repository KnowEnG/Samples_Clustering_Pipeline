import os
import unittest
from unittest import TestCase
import numpy as np

import knpackage.toolbox as kn
import sample_clustering_toolbox_research_module as tstdata
import sample_clustering_toolbox as sctbx

class TestSave_consensus_clustering(TestCase):
    def setUp(self):
        self.run_parameters = tstdata.get_test_paramters_dictionary()
        self.run_parameters["run_directory"] = '.'
        self.run_parameters["tmp_directory"] = kn.create_dir(self.run_parameters["run_directory"], 'tmp_cc_net_nmf')
        self.run_parameters["results_directory"] = self.run_parameters["tmp_directory"]
        self.run_parameters['display_clusters'] = 0
        self.run_parameters["use_now_name"] = 0

    def tearDown(self):
        kn.remove_dir(self.run_parameters["tmp_directory"])
        del self.run_parameters

    def test_save_consensus_clustering(self):
        n_clusters = int(self.run_parameters['number_of_clusters'])
        spreadsheet_consensus_mat = tstdata.get_square_3_cluster_spreadsheet(n_clusters)
        n_rows = spreadsheet_consensus_mat.shape[0]
        sample_names = np.random.permutation(n_rows)
        labels = np.random.randint(0, n_clusters, n_rows)

        sctbx.save_consensus_clustering(spreadsheet_consensus_mat, sample_names, labels, self.run_parameters)
        labels_data_file = self.run_parameters['consensus_clustering_file']

        self.assertTrue(os.path.exists(labels_data_file), msg='nmf - labels_data.tsv:  dne')

if __name__ == '__main__':
    unittest.main()