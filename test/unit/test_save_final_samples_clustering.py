import os
import unittest
from unittest import TestCase
import numpy as np

import knpackage.toolbox as kn
import sample_clustering_toolbox_research_module as tstdata
import sample_clustering_toolbox as sctbx

class TestSave_final_samples_clustering(TestCase):
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

    def test_save_final_samples_clustering(self):
        n_clusters = int(self.run_parameters['number_of_clusters'])
        n_rows = 100
        sample_names = np.random.permutation(n_rows)
        labels = np.random.randint(0, n_clusters, n_rows)
        sctbx.save_final_samples_clustering(sample_names, labels, self.run_parameters)
        labels_data_file = self.run_parameters['cluster_labels_file']

        self.assertTrue(os.path.exists(labels_data_file), msg='nmf - labels_data.tsv:  dne')

if __name__ == '__main__':
    unittest.main()

"""
def save_final_samples_clustering(sample_names, labels, run_parameters):
    #"" wtite .tsv file that assings a cluster number label to the sample_names.

    Args:
        sample_names: (unique) data identifiers.
        labels: cluster number assignments.
        run_parameters: write path (run_parameters["results_directory"]).
    #""
    file_name = os.path.join(run_parameters["results_directory"], kn.create_timestamped_filename('labels_data', 'tsv'))
    df_tmp = kn.create_df_with_sample_labels(sample_names, labels)
    df_tmp.to_csv(file_name, sep='\t', header=None)
    run_parameters['cluster_labels_file'] = file_name

"""