# -*- coding: utf-8 -*-
"""
Demonstrate two ways to create and run a unittest.TestSuite

Created on Tue May 31 13:50:35 2016
Python Script for Unit Test of KnowEnG.py modules
KnowEnG_UnitTestsRunner.py

commandLine>> run KnowEnG_UnitTestsRunner.py

lanier4@illinois.edu
"""
sample_clustering_unit_test_directory = '/Users/lanier4/PycharmProjects/Samples_Clustering_Pipeline/test/unit'
import unittest
import sys
sys.path.extend(sample_clustering_unit_test_directory)
import test_sample_clustering_toolbox as tkeg

#                               method B
mySuit = tkeg.suite()
runner = unittest.TextTestRunner()
myResult2 = runner.run(mySuit)