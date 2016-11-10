import os
import unittest
from unittest import TestCase

import sample_clustering_toolbox_research_module as tstdata
import sample_clustering_toolbox as sctbx

class TestRun_nmf(TestCase):
    def test_run_nmf(self):
        run_parameters = tstdata.get_test_paramters_dictionary()

        run_parameters["use_now_name"] = 0

        sctbx.run_nmf(run_parameters)
        file_name = run_parameters['cluster_labels_file']
        self.assertTrue(os.path.exists(file_name), msg='nmf - labels_data.tsv:  dne')

if __name__ == '__main__':
    unittest.main()