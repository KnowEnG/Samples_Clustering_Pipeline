# -*- coding: utf-8 -*-
"""     lanier4@illinois.edu    """

import unittest
import numpy as np
import numpy.linalg as lpac
import scipy.sparse as spar

# file version to test
import knpackage.toolbox as keg

class toolbox_test(unittest.TestCase):
    
    def testQN(self):
        a = np.array([[7.0, 5.0],[3.0, 1.0],[1.0,7.0]])
        aQN = np.array([[7.0, 4.0],[4.0,1.0],[1.0,7.0]])
        qn1 = keg.get_quantile_norm_matrix(a)
        
        self.assertEqual(sum(sum(qn1 != aQN)), 0, 'Quantile Norm 1 Not Equal')
        
    def testQN_T(self):
        a = np.array([[7.0, 5.0],[3.0, 1.0],[1.0,7.0]])
        aQN = np.array([[7.0, 4.0],[4.0,1.0],[1.0,7.0]])
        qn1 = keg.get_quantile_norm_matrix(a.T)
        qn1 = qn1.T
        
        self.assertEqual(sum(sum(qn1 != aQN)), 0, 'Quantile Norm 1 Not Equal')
        
    def testRWR(self):

        alpha = 1.0
        maxIt = 100
        tol = 1e-8

        F0 = np.array([1.0, 0.0])
        A = spar.csr_matrix(( (np.eye(2) + np.ones(2)) / 3) )
        A = spar.csr_matrix(( (np.eye(2) + np.ones(2)) / 3) )
        
        F_exact = np.array([0.5, 0.5])
        F_calculated = keg.smooth_matrix_with_rwr(F0, A, alpha, maxIt, tol)
        self.assertLessEqual( lpac.norm(F_exact - F_calculated), tol)
        
    def testRWR_2(self):

        alpha = 0.5
        maxIt = 100
        tol = 1e-8

        F0 = np.array([1.0, 0.0])
        A = spar.csr_matrix(((np.eye(2) + np.ones(2)) / 3))
        #A = spar.csr_matrix(A)
        
        F_exact = np.array([0.8, 0.2])
        F_calculated = keg.smooth_matrix_with_rwr(F0, A, alpha, maxIt, tol)
        self.assertLessEqual( lpac.norm(F_exact - F_calculated), tol)     
        
        

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