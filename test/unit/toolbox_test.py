# -*- coding: utf-8 -*-
"""     lanier4@illinois.edu    """

import unittest
import os
import numpy as np
#import numpy.linalg as LA
import scipy.sparse as spar
import scipy.stats as stats

# file version to test
import knpackage.toolbox as kn

class toolbox_test(unittest.TestCase):
    
    def get_run_parameters(self):
        run_parameters = {'test_directory':'/Users/lanier4/BigDataTank/unit_test_development',
                            'k':3,'number_of_iteriations_in_rwr':100,
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
    
    def test_get_quantile_norm_matrix(self):
        a = np.array([[7.0, 5.0],[3.0, 1.0],[1.0,7.0]])
        aQN = np.array([[7.0, 4.0],[4.0,1.0],[1.0,7.0]])
        qn1 = kn.get_quantile_norm_matrix(a)
        
        self.assertEqual(sum(sum(qn1 != aQN)), 0, 'Quantile Norm 1 Not Equal')
        
    def test_smooth_matrix_with_rwr(self):
        """ Assert that a test matrix will converge to the precomputed answer in
            the predicted number of steps (iterations). Depends on run_parameters
            and the values set herein.
        """
        EXPECTED_STEPS = 25
        run_parameters = self.get_run_parameters()
        run_parameters['restart_probability'] = 1
        run_parameters['restart_tolerance'] = 1e-12
        F0 = np.eye(2)
        A = spar.csr_matrix(( (np.eye(2) + np.ones(2)) / 3) )
        # A = (np.eye(2) + np.ones((2,2))) / 3
        
        F_exact = np.ones((2, 2)) * 0.5
        F_calculated, steps = kn.smooth_matrix_with_rwr(F0, A, run_parameters)
        self.assertEqual(steps, EXPECTED_STEPS)
        
        T = (np.abs(F_exact - F_calculated))
        
        self.assertAlmostEqual(T.sum(), 0)
        
    def test_smooth_matrix_with_rwr_non_sparse(self):
        """ Assert that a test matrix will converge to the precomputed answer in
            the predicted number of steps (iterations). Depends on run_parameters
            and the values set herein.
        """
        EXPECTED_STEPS = 25
        run_parameters = self.get_run_parameters()
        run_parameters['restart_probability'] = 1
        run_parameters['restart_tolerance'] = 1e-12
        F0 = np.eye(2)
        #A = spar.csr_matrix(( (np.eye(2) + np.ones(2)) / 3) )
        A = (np.eye(2) + np.ones((2,2))) / 3
        
        F_exact = np.ones((2, 2)) * 0.5
        F_calculated, steps = kn.smooth_matrix_with_rwr(F0, A, run_parameters)
        self.assertEqual(steps, EXPECTED_STEPS)
        
        T = (np.abs(F_exact - F_calculated))
        
        self.assertAlmostEqual(T.sum(), 0)
        
    def test_smooth_matrix_with_rwr_single_vector(self):
        """ Assert that a test matrix will converge to the precomputed answer in
            the predicted number of steps (iterations). Depends on run_parameters
            and the values set herein.
        """
        EXPECTED_STEPS = 25
        run_parameters = self.get_run_parameters()
        run_parameters['restart_probability'] = 1
        run_parameters['restart_tolerance'] = 1e-12
        F0 = np.array([1.0, 0.0])
        #A = spar.csr_matrix(( (np.eye(2) + np.ones(2)) / 3) )
        A = (np.eye(2) + np.ones((2,2))) / 3
        
        F_exact = np.array([0.5, 0.5])
        F_calculated, steps = kn.smooth_matrix_with_rwr(F0, A, run_parameters)
        self.assertEqual(steps, EXPECTED_STEPS, msg='minor difference')
        T = (np.abs(F_exact - F_calculated))
        self.assertAlmostEqual(T.sum(), 0)
        
    def test_normalize_sparse_mat_by_diagonal(self):
        """ assert that a test matrix will be "normalized" s.t. the sum of the rows
            or columns will nearly equal one
        """
        A = np.random.rand(500,500)
        B = kn.normalize_sparse_mat_by_diagonal(spar.csr_matrix(A))
        B = B.todense()
        B2 = B**2
        B2 = np.sqrt(B2.sum())
        geo_mean = float(stats.gmean(B.sum(axis=1)))
        self.assertAlmostEqual(geo_mean, 1, delta=0.1)
        geo_mean = float(stats.gmean(B.T.sum(axis=1)))
        self.assertAlmostEqual(geo_mean, 1, delta=0.1)
        
    def test_form_network_laplacian_matrix(self):
        """ assert that the laplacian matrix returned sums to zero in both rows
            and columns
        """
        THRESHOLD = 0.8
        A  = np.random.rand(10,10)
        A[A < THRESHOLD] = 0
        A = A + A.T
        Ld, Lk = kn.form_network_laplacian_matrix(A)
        L = Ld - Lk
        L = L.todense()
        L0 = L.sum(axis=0)
        self.assertFalse(L0.any(), msg='Laplacian row sum not equal 0')
        L1 = L.sum(axis=1)
        self.assertFalse(L1.any(), msg='Laplacian col sum not equal 0')
        Ld, Lk = kn.form_network_laplacian_matrix(spar.csr_matrix(A))
        L = Ld - Lk
        L = L.todense()
        L0 = L.sum(axis=0)
        self.assertFalse(L0.any(), msg='Laplacian row sum not equal 0')
        L1 = L.sum(axis=1)
        self.assertFalse(L1.any(), msg='Laplacian col sum not equal 0')
        
    def test_sample_a_matrix(self):
        """ assert that the random sample is of the propper size, the
            permutation points to the correct columns and that the number of 
            rows set to zero is correct.
        """
        n_test_rows = 11
        n_test_cols = 5
        pct_smpl = 0.6
        n_zero_rows = int(np.round(n_test_rows * (1 - pct_smpl)))
        n_smpl_cols = int(np.round(n_test_cols * pct_smpl))
        epsilon_sum = max(n_test_rows, n_test_cols) * 1e-15
        A = np.random.rand(n_test_rows, n_test_cols) + epsilon_sum
        B, P = kn.sample_a_matrix(A, pct_smpl)
        self.assertEqual(B.shape[1], P.size, msg='permutation size not equal columns')
        self.assertEqual(P.size, n_smpl_cols, msg='number of sample columns exception')
        perm_err_sum = 0
        n_zero_err_sum = 0
        B_col = 0
        for A_col in P:
            n_zeros = (np.int_(B[:, B_col] == 0)).sum()
            if n_zeros != n_zero_rows:
                n_zero_err_sum += 1
            C = A[:, A_col] - B[:, B_col]
            C[B[:, B_col] == 0] = 0
            B_col += 1
            if C.sum() > epsilon_sum:
                perm_err_sum += 1
        
        self.assertEqual(n_zero_err_sum, 0, msg='number of zero columns exception')
        self.assertEqual(perm_err_sum, 0, msg='permutation index exception')

    def test_create_dir_AND_remove_dir(self):
        """ assert that the functions work togeather to create and remove a directory
            even when files have been added
        """
        dir_name = 'tmp_test'
        run_parameters = self.get_run_parameters()
        dir_path = run_parameters['test_directory']
        ndr = kn.create_dir(dir_path, dir_name)
        self.assertTrue(os.path.exists(ndr), msg='create_dir function exception')
        A = np.random.rand(10,10)
        time_stamp = '123456789'
        a_name = os.path.join(ndr, 'temp_test' + time_stamp)
        A.dump(a_name)
        A_back = np.load(a_name)
        A_diff = A - A_back
        A_diff = A_diff.sum()
        self.assertEqual(A_diff, 0, msg='write / read directory exception')
        kn.remove_dir(ndr)
        self.assertFalse(os.path.exists(ndr), msg='remove_dir function exception')
        
    def test_get_timestamp(self):
        """ assert that the default size of the timestamp string is 16 chars and
            that sequential calls produce differnt results
        """
        n_default_chars = 16
        stamp_time = 1e6
        tstr = kn.get_timestamp(stamp_time)
        tstr2 = kn.get_timestamp()
        
        self.assertEqual(len(tstr), n_default_chars, msg='string return size unexpected')
        self.assertNotEqual(tstr, tstr2)
        
        
    """
def get_timestamp(stamp_units=1e6):
    get a time stamp string - current time as integer string.

    Args:
        stamp_units: inverse of time resolution 1e6 returns microseconds.

    Returns:
        timestamp_string: a string of integer digits.

    timestamp_string = np.str_(int(time.time() * np.maximum(stamp_units, 1)))

    return timestamp_string

def create_timestamped_filename(name_base='t', stamp_units=1e6):
    append a filename with a timestamp string.

    Args:
        name_base: the file name - a prefix to the time stamp string.
        stamp_units: time resolution; 1e6 for microseconds, 1e3 milliseconds.

    Returns:
        time_stamped_file_name: name_base_123456 (some long number)

    time_stamped_file_name = name_base + '_' + get_timestamp(stamp_units)

    return time_stamped_file_name

def append_run_parameters_dict(run_parameters, key_name, value_str):
    add a key-value pair to the run parameters dictionary.

    Args:
        run_parameters: dictionary to append.
        key_name: key name to add or overwrite.
        value_str: value to insert in run_parameters[key_name].

    Returns:
        run_parameters: dictionary with new (or overwritten) key value pair.

    run_parameters[key_name] = value_str

    return run_parameters    
    
def update_h_coordinate_matrix(w_matrix, x_matrix):
    nonnegative right factor matrix for perform_net_nmf function s.t. X ~ W.H.

    Args:
        w_matrix: the positive left factor (W) of the perform_net_nmf function.
        x_matrix: the postive matrix (X) to be decomposed.

    Returns:
        h_matrix: nonnegative right factor (H) matrix.

    wtw = np.dot(w_matrix.T, w_matrix)
    number_of_clusters = wtw.shape[0]
    wtx = np.dot(w_matrix.T, x_matrix)
    colix = np.arange(0, x_matrix.shape[1])
    rowix = np.arange(0, w_matrix.shape[1])
    h_matrix = np.dot(LA.pinv(wtw), wtx)
    h_pos = h_matrix > 0
    h_matrix[~h_pos] = 0
    col_log_arr = sum(h_pos == 0) > 0
    col_list = colix[col_log_arr]
    for cluster in range(0, number_of_clusters):
        if col_list.size > 0:
            w_ette = wtx[:, col_list]
            m_rows = w_ette.shape[0]
            n_cols = w_ette.shape[1]
            mcode_uniq_col_ix = np.arange(0, n_cols)
            h_ette = np.zeros((m_rows, n_cols))
            h_pos_ette = h_pos[:, col_list]
            mcoding = np.dot(2**(np.arange(0, m_rows)), np.int_(h_pos_ette))
            mcode_uniq = np.unique(mcoding)
            for u_n in mcode_uniq:
                ixidx = mcoding == u_n
                c_pat = mcode_uniq_col_ix[ixidx]
                if c_pat.size > 0:
                    r_pat = rowix[h_pos_ette[:, c_pat[0]]]
                    atmp = wtw[r_pat[:, None], r_pat]
                    btmp = w_ette[r_pat[:, None], c_pat]
                    atmptatmp = np.dot(atmp.T, atmp)
                    atmptatmp = LA.pinv(atmptatmp)
                    atmptbtmp = np.dot(atmp.T, btmp)
                    h_ette[r_pat[:, None], c_pat] = np.dot(atmptatmp, atmptbtmp)
                    h_matrix[:, col_list] = h_ette
            h_pos = h_matrix > 0
            h_matrix[~h_pos] = 0
            col_log_arr = sum(h_pos == 0) > 0
            col_list = colix[col_log_arr]
        else:
            break

    return h_matrix
    
def perform_net_nmf(x_matrix, lap_val, lap_dag, run_parameters):
    perform network based nonnegative matrix factorization, minimize:
        ||X-WH|| + lambda.tr(W'.L.W), with W, H positive.

    Args:
        x_matrix: the postive matrix (X) to be decomposed into W.H
        lap_val: the laplacian matrix
        lap_dag: the diagonal of the laplacian matrix
        run_parameters: parameters dictionary with keys: "k", "lambda", "it_max",
            "h_clust_eq_limit", "obj_fcn_chk_freq".

    Returns:
        h_matrix: nonnegative right factor (H) matrix.

    k = int(run_parameters["k"])
    lmbda = float(run_parameters["lmbda"])
    epsilon = 1e-15
    w_matrix = np.random.rand(x_matrix.shape[0], k)
    w_matrix = maximum(w_matrix / maximum(sum(w_matrix), epsilon), epsilon)
    h_matrix = np.random.rand(k, x_matrix.shape[1])
    h_clust_eq = np.argmax(h_matrix, 0)
    h_eq_count = 0
    for itr in range(0, int(run_parameters["it_max"])):
        if np.mod(itr, int(run_parameters["obj_fcn_chk_freq"])) == 0:
            h_clusters = np.argmax(h_matrix, 0)
            if (itr > 0) & (sum(h_clust_eq != h_clusters) == 0):
                h_eq_count = h_eq_count + int(run_parameters["obj_fcn_chk_freq"])
            else:
                h_eq_count = 0
            h_clust_eq = h_clusters
            if h_eq_count >= float(run_parameters["h_clust_eq_limit"]):
                break
        numerator = maximum(np.dot(x_matrix, h_matrix.T) + lmbda * lap_val.dot(w_matrix), epsilon)
        denomerator = maximum(np.dot(w_matrix, np.dot(h_matrix, h_matrix.T))
                              + lmbda * lap_dag.dot(w_matrix), epsilon)
        w_matrix = w_matrix * (numerator / denomerator)
        w_matrix = maximum(w_matrix / maximum(sum(w_matrix), epsilon), epsilon)
        h_matrix = update_h_coordinate_matrix(w_matrix, x_matrix)

    return h_matrix

def perform_nmf(x_matrix, run_parameters):
    nonnegative matrix factorization, minimize the diffence between X and W dot H
        with positive factor matrices W, and H.

    Args:
        x_matrix: the postive matrix (X) to be decomposed into W dot H.
        run_parameters: parameters dictionary with keys "k", "it_max",
            "cluster_min_repeats", "obj_fcn_chk_freq".

    Returns:
        h_matrix: nonnegative right factor matrix (H).

    k = int(run_parameters["k"])
    obj_fcn_chk_freq = int(run_parameters["obj_fcn_chk_freq"])
    h_clust_eq_limit = float(run_parameters["h_clust_eq_limit"])
    epsilon = 1e-15
    w_matrix = np.random.rand(x_matrix.shape[0], k)
    w_matrix = maximum(w_matrix / maximum(sum(w_matrix), epsilon), epsilon)
    h_matrix = np.random.rand(k, x_matrix.shape[1])
    h_clust_eq = np.argmax(h_matrix, 0)
    h_eq_count = 0
    for itr in range(0, int(run_parameters["it_max"])):
        if np.mod(itr, obj_fcn_chk_freq) == 0:
            h_clusters = np.argmax(h_matrix, 0)
            if (itr > 0) & (sum(h_clust_eq != h_clusters) == 0):
                h_eq_count = h_eq_count + obj_fcn_chk_freq
            else:
                h_eq_count = 0
            h_clust_eq = h_clusters
            if h_eq_count >= h_clust_eq_limit:
                break
        numerator = maximum(np.dot(x_matrix, h_matrix.T), epsilon)
        denomerator = maximum(np.dot(w_matrix, np.dot(h_matrix, h_matrix.T)), epsilon)
        w_matrix = w_matrix * (numerator / denomerator)
        w_matrix = maximum(w_matrix / maximum(sum(w_matrix), epsilon), epsilon)
        h_matrix = update_h_coordinate_matrix(w_matrix, x_matrix)

    return h_matrix

def update_linkage_matrix(encode_mat, sample_perm, linkage_matrix):
    update the connectivity matrix by summing the un-permuted linkages.
    encode_mat: (permuted) nonnegative right factor matrix (H) - encoded linkage.
    Args:
        encode_mat: encoding of linkage either as an h_matrix or argmax(h_matrix)

        sample_perm: the sample permutaion of the h_matrix.
        linkage_matrix: connectivity matrix.

    Returns:
        linkage_matrix: connectivity matrix summed with the de-permuted linkage.

    if encode_mat.ndim == 1:
        num_clusters = max(encode_mat) + 1
        cluster_id = encode_mat
    else:
        num_clusters = encode_mat.shape[0]
        cluster_id = np.argmax(encode_mat, 0)

    for cluster in range(0, num_clusters):
        slice_id = sample_perm[cluster_id == cluster]
        linkage_matrix[slice_id[:, None], slice_id] += 1

    return linkage_matrix

def update_indicator_matrix(sample_perm, indicator_matrix):
    update the indicator matrix by summing the un-permutation.

    Args:
        sample_perm: permutaion of the sample (h_matrix).
        indicator_matrix: indicator matrix.

    Returns:
        indicator_matrix: indicator matrix incremented at sample_perm locations.

    indicator_matrix[sample_perm[:, None], sample_perm] += 1

    return indicator_matrix

def perform_kmeans(consensus_matrix, k=3):
    determine cluster assignments for consensus matrix using K-means.

    Args:
        consensus_matrix: connectivity / indicator matrix.
        k: clusters estimate.

    Returns:
        lablels: ordered cluster assignments for consensus_matrix (samples).

    cluster_handle = KMeans(k, random_state=10)
    labels = cluster_handle.fit_predict(consensus_matrix)

    return labels
    """
        

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
mySuit = tkn.suite()
runner = unittest.TextTestRunner()
myResult = runner.run(mySuit)

OR
mySuit2 = unittest.TestLoader().loadTestsFromTestCase(TestKEGmodule)

'''    