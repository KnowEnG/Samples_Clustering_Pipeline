import os
import unittest
from unittest import TestCase

import sample_clustering_toolbox_research_module as tstdata
import sample_clustering_toolbox as sctbx

class TestRun_cc_net_nmf(TestCase):

    def test_run_cc_net_nmf(self):
        run_parameters = tstdata.get_test_paramters_dictionary()
        run_parameters['method'] = run_parameters['method_4']

        run_parameters['display_clusters'] = 0
        run_parameters["use_now_name"] = 0
        run_parameters["processing_method"] = 'parl_loc' 

        sctbx.run_cc_net_nmf(run_parameters)
        file_name = run_parameters['cluster_labels_file']
        self.assertTrue(os.path.exists(file_name), msg='nmf - labels_data.tsv:  dne')

if __name__ == '__main__':
    unittest.main()
