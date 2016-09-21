#import os
import unittest
from unittest import TestCase
import numpy as np

#import knpackage.toolbox as kn
import sample_clustering_toolbox as sctbx

class TestPerform_pearson_correlation(TestCase):

    def test_perform_pearson_correlation(self):
        spreadsheet = np.array([[0.1, -0.2], [-0.3, 0.6], [0.4, -0.1]])
        drug_response = np.array([0.5, 0.7])
        expected_pc = np.array([-1.0, 1.0, -1.0])
        pc_array = sctbx.perform_pearson_correlation(spreadsheet, drug_response)
        for n in range(0, expected_pc.size):
            self.assertEqual(expected_pc[n], pc_array[n], msg='pearson coefficent unexpected')

if __name__ == '__main__':
    unittest.main()