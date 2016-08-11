# -*- coding: utf-8 -*-
"""     lanier4@illinois.edu    """

import unittest
import numpy as np
#import numpy.linalg as LA
import scipy.sparse as spar

# file version to test
import knpackage.toolbox as keg

class toolbox_test(unittest.TestCase):
    
    def get_run_parameters(self):
        run_parameters = {'k':3,'number_of_iteriations_in_rwr':100,
                          'obj_fcn_chk_freq':50,
                          'it_max':10000,
                          'h_clust_eq_limit':100,
                          'restart_tolerance':0.0001,
                          'lmbda':1400,
                          'percent_sample':0.8,
                          'number_of_bootstraps':3,
                          'display_clusters':1,
                          'restart_probability':0.7,
                          'verbose':1,
                          'use_now_name':1000000}

        return run_parameters
    
    def testQN(self):
        a = np.array([[7.0, 5.0],[3.0, 1.0],[1.0,7.0]])
        aQN = np.array([[7.0, 4.0],[4.0,1.0],[1.0,7.0]])
        qn1 = keg.get_quantile_norm_matrix(a)
        
        self.assertEqual(sum(sum(qn1 != aQN)), 0, 'Quantile Norm 1 Not Equal')
        
    def testRWR(self):
        run_parameters = self.get_run_parameters()
        run_parameters['restart_probability'] = 1
        run_parameters['restart_tolerance'] = 1e-12
        F0 = np.eye(2)
        A = spar.csr_matrix(( (np.eye(2) + np.ones(2)) / 3) )
        # A = (np.eye(2) + np.ones((2,2))) / 3
        
        F_exact = np.ones((2, 2)) * 0.5
        F_calculated, steps = keg.smooth_matrix_with_rwr(F0, A, run_parameters)
        T = (np.abs(F_exact - F_calculated))
        
        self.assertAlmostEqual(T.sum(), 0)
        

def suite():
    test_suite = unittest.TestSuite()
    test_suite.addTest(unittest.makeSuite(toolbox_test))
    
    return test_suite

'''# Next two lines for using this file w/o test Suite   << NOT recommended
#if __name__=='__main__':
#    unittest.main()

                                        >> Preferred Method for using unit test
import unittest
import TestKEGmodule as tkeg
mySuit = tkeg.suite()
runner = unittest.TextTestRunner()
myResult = runner.run(mySuit)

OR
mySuit2 = unittest.TestLoader().loadTestsFromTestCase(TestKEGmodule)

'''    