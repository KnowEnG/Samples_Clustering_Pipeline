import os
import unittest
from unittest import TestCase

import test.unit.sample_clustering_toolbox_research_module as tstdata
import src.sample_clustering_toolbox as sctbx

class TestRun_net_nmf(TestCase):
    def test_run_net_nmf(self):
        run_parameters = tstdata.get_test_paramters_dictionary()
        run_parameters['method'] = run_parameters['method_3']

        run_parameters['display_clusters'] = 0
        run_parameters["use_now_name"] = 0

        sctbx.run_net_nmf(run_parameters)
        file_name = run_parameters['cluster_labels_file']
        self.assertTrue(os.path.exists(file_name), msg='nmf - labels_data.tsv:  dne')

if __name__ == '__main__':
    unittest.main()
