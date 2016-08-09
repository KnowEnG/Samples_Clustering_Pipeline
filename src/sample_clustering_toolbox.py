# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 16:08:25 2016

@author: del
@author: The Gene Sets Characterization dev team

"""
import os
import time
#import argparse
import numpy as np
#import numpy.linalg as LA
#from numpy import maximum

import pandas as pd
#import scipy.sparse as spar
#from scipy import stats
#from sklearn.preprocessing import normalize
#from sklearn.cluster import KMeans
import yaml
import matplotlib.pyplot as plt

import knpackage.toolbox as keg

def run_nmf(run_parameters):
    """ wrapper: call sequence to perform non-negative matrix factorization and write results.

    Args:
        run_parameters: parameter set dictionary.
    """
    spreadsheet_df = keg.get_spreadsheet_df(run_parameters)
    spreadsheet_mat = spreadsheet_df.as_matrix()
    spreadsheet_mat = keg.get_quantile_norm_matrix(spreadsheet_mat)

    h_mat = keg.perform_nmf(spreadsheet_mat, run_parameters)

    linkage_matrix = np.zeros((spreadsheet_mat.shape[1], spreadsheet_mat.shape[1]))
    sample_perm = np.arange(0, spreadsheet_mat.shape[1])
    linkage_matrix = keg.update_linkage_matrix(h_mat, sample_perm, linkage_matrix)
    labels = keg.perform_kmeans(linkage_matrix, int(run_parameters['k']))

    sample_names = spreadsheet_df.columns
    save_clusters(sample_names, labels, run_parameters)

    if int(run_parameters['display_clusters']) != 0:
        con_mat_image = form_consensus_matrix_graphic(linkage_matrix, int(run_parameters['k']))
        display_clusters(con_mat_image)

    return

def run_cc_nmf(run_parameters):
    """ wrapper: call sequence to perform non-negative matrix factorization with
        consensus clustering and write results.

    Args:
        run_parameters: parameter set dictionary.
    """
    spreadsheet_df = keg.get_spreadsheet_df(run_parameters)
    spreadsheet_mat = spreadsheet_df.as_matrix()
    spreadsheet_mat = keg.get_quantile_norm_matrix(spreadsheet_mat)

    find_and_save_nmf_clusters(spreadsheet_mat, run_parameters)

    linkage_matrix = np.zeros((spreadsheet_mat.shape[1], spreadsheet_mat.shape[1]))
    indicator_matrix = linkage_matrix.copy()
    consensus_matrix = form_consensus_matrix(run_parameters, linkage_matrix, indicator_matrix)
    labels = keg.perform_kmeans(consensus_matrix, int(run_parameters['k']))

    sample_names = spreadsheet_df.columns
    save_consensus_cluster_result(consensus_matrix, sample_names, labels, run_parameters)

    if int(run_parameters['display_clusters']) != 0:
        display_clusters(form_consensus_matrix_graphic(consensus_matrix, int(run_parameters['k'])))

    return

def run_net_nmf(run_parameters):
    """ wrapper: call sequence to perform network based stratification and write results.

    Args:
        run_parameters: parameter set dictionary.
    """
    spreadsheet_df = keg.get_spreadsheet_df(run_parameters)
    network_df = keg.get_network_df(run_parameters['network_file_name'])

    node_1_names, node_2_names = keg.extract_network_node_names(network_df)
    unique_gene_names = keg.find_unique_node_names(node_1_names, node_2_names)
    genes_lookup_table = keg.create_node_names_dict(unique_gene_names)

    network_df = keg.map_node_names_to_index(network_df, genes_lookup_table, 'node_1')
    network_df = keg.map_node_names_to_index(network_df, genes_lookup_table, 'node_2')

    network_df = keg.symmetrize_df(network_df)
    network_mat = keg.convert_network_df_to_sparse(network_df, len(unique_gene_names), len(unique_gene_names))

    network_mat = keg.normalized_matrix(network_mat)
    lap_diag, lap_pos = keg.form_network_laplacian_matrix(network_mat)

    spreadsheet_df = keg.update_spreadsheet_df(spreadsheet_df, unique_gene_names)
    spreadsheet_mat = spreadsheet_df.as_matrix()
    sample_names = spreadsheet_df.columns

    sample_smooth, iterations = keg.smooth_matrix_with_rwr(
        spreadsheet_mat, network_mat, run_parameters)
    sample_quantile_norm = keg.get_quantile_norm_matrix(sample_smooth)
    h_mat = keg.perform_net_nmf(sample_quantile_norm, lap_pos, lap_diag, run_parameters)

    linkage_matrix = np.zeros((spreadsheet_mat.shape[1], spreadsheet_mat.shape[1]))
    sample_perm = np.arange(0, spreadsheet_mat.shape[1])
    linkage_matrix = keg.update_linkage_matrix(h_mat, sample_perm, linkage_matrix)
    labels = keg.perform_kmeans(linkage_matrix, int(run_parameters["k"]))

    save_clusters(sample_names, labels, run_parameters)

    if int(run_parameters['display_clusters']) != 0:
        display_clusters(form_consensus_matrix_graphic(linkage_matrix, int(run_parameters['k'])))

    return

def run_cc_net_nmf(run_parameters):
    """ wrapper: call sequence to perform network based stratification with consensus clustering
        and write results.

    Args:
        run_parameters: parameter set dictionary.
    """
    spreadsheet_df = keg.get_spreadsheet_df(run_parameters)
    network_df = keg.get_network_df(run_parameters['network_file_name'])

    node_1_names, node_2_names = keg.extract_network_node_names(network_df)
    unique_gene_names = keg.find_unique_node_names(node_1_names, node_2_names)
    genes_lookup_table = keg.create_node_names_dict(unique_gene_names)

    network_df = keg.map_node_names_to_index(network_df, genes_lookup_table, 'node_1')
    network_df = keg.map_node_names_to_index(network_df, genes_lookup_table, 'node_2')

    network_df = keg.symmetrize_df(network_df)
    #network_mat = convert_df_to_sparse(network_df, len(unique_gene_names))
    network_mat = keg.convert_network_df_to_sparse(network_df, len(unique_gene_names), len(unique_gene_names))

    network_mat = keg.normalized_matrix(network_mat)
    lap_diag, lap_pos = keg.form_network_laplacian_matrix(network_mat)

    spreadsheet_df = keg.update_spreadsheet_df(spreadsheet_df, unique_gene_names)
    spreadsheet_mat = spreadsheet_df.as_matrix()
    sample_names = spreadsheet_df.columns

    find_and_save_net_nmf_clusters(network_mat, spreadsheet_mat, lap_diag, lap_pos, run_parameters)

    linkage_matrix = np.zeros((spreadsheet_mat.shape[1], spreadsheet_mat.shape[1]))
    indicator_matrix = linkage_matrix.copy()
    consensus_matrix = form_consensus_matrix(run_parameters, linkage_matrix, indicator_matrix)
    labels = keg.perform_kmeans(consensus_matrix, int(run_parameters['k']))

    save_consensus_cluster_result(consensus_matrix, sample_names, labels, run_parameters)

    if int(run_parameters['display_clusters']) != 0:
        display_clusters(form_consensus_matrix_graphic(consensus_matrix, int(run_parameters['k'])))

    return
    
def find_and_save_net_nmf_clusters(network_mat, spreadsheet_mat, lap_dag, lap_val, run_parameters):
    """ central loop: compute components for the consensus matrix from the input
        network and spreadsheet matrices and save them to temp files.

    Args:
        network_mat: genes x genes symmetric matrix.
        spreadsheet_mat: genes x samples matrix.
        lap_dag, lap_val: laplacian matrix components; L = lap_dag - lap_val.
        run_parameters: dictionay of run-time parameters.
    """
    for sample in range(0, int(run_parameters["number_of_bootstraps"])):
        sample_random, sample_permutation = keg.sample_a_matrix(
            spreadsheet_mat, np.float64(run_parameters["percent_sample"]))
        sample_smooth, iterations = \
        keg.smooth_matrix_with_rwr(sample_random, network_mat, run_parameters)

        if int(run_parameters['verbose']) != 0:
            print("{} of {}: iterations = {}".format(
                sample + 1, run_parameters["number_of_bootstraps"], iterations))

        sample_quantile_norm = keg.get_quantile_norm_matrix(sample_smooth)
        h_mat = keg.perform_net_nmf(sample_quantile_norm, lap_val, lap_dag, run_parameters)

        save_temporary_cluster(h_mat, sample_permutation, run_parameters, sample)

    return

def find_and_save_nmf_clusters(spreadsheet_mat, run_parameters):
    """ central loop: compute components for the consensus matrix by
        non-negative matrix factorization.

    Args:
        spreadsheet_mat: genes x samples matrix.
        run_parameters: dictionay of run-time parameters.
    """
    for sample in range(0, int(run_parameters["number_of_bootstraps"])):
        sample_random, sample_permutation = keg.sample_a_matrix(
            spreadsheet_mat, np.float64(run_parameters["percent_sample"]))

        h_mat = keg.perform_nmf(sample_random, run_parameters)
        save_temporary_cluster(h_mat, sample_permutation, run_parameters, sample)

        if int(run_parameters['verbose']) != 0:
            print('nmf {} of {}'.format(
                sample + 1, run_parameters["number_of_bootstraps"]))

    return
    
def form_consensus_matrix(run_parameters, linkage_matrix, indicator_matrix):
    """ compute the consensus matrix from the indicator and linkage matrix inputs
        formed by the bootstrap "temp_*" files.

    Args:
        run_parameters: parameter set dictionary with "tmp_directory" key.
        linkage_matrix: linkage matrix from initialization or previous call.
        indicator_matrix: indicator matrix from initialization or previous call.

    Returns:
        consensus_matrix: (sum of linkage matrices) / (sum of indicator matrices).
    """
    indicator_matrix = form_indicator_matrix(run_parameters, indicator_matrix)
    linkage_matrix = form_linkage_matrix(run_parameters, linkage_matrix)
    consensus_matrix = linkage_matrix / np.maximum(indicator_matrix, 1)

    return consensus_matrix

def form_indicator_matrix(run_parameters, indicator_matrix):
    """ read bootstrap temp_p* files and compute the indicator_matrix.

    Args:
        run_parameters: parameter set dictionary.
        indicator_matrix: indicator matrix from initialization or previous call.

    Returns:
        indicator_matrix: input summed with "temp_p*" files in run_parameters["tmp_directory"].
    """
    tmp_dir = run_parameters["tmp_directory"]
    dir_list = os.listdir(tmp_dir)
    for tmp_f in dir_list:
        if tmp_f[0:6] == 'temp_p':
            pname = os.path.join(tmp_dir, tmp_f)
            sample_permutation = np.load(pname)
            indicator_matrix = keg.update_indicator_matrix(sample_permutation, indicator_matrix)

    return indicator_matrix

def form_linkage_matrix(run_parameters, linkage_matrix):
    """ read bootstrap temp_h* and temp_p* files, compute and add the linkage_matrix.

    Args:
        run_parameters: parameter set dictionary.
        linkage_matrix: connectivity matrix from initialization or previous call.

    Returns:
        linkage_matrix: summed with "temp_h*" files in run_parameters["tmp_directory"].
    """
    tmp_dir = run_parameters["tmp_directory"]
    dir_list = os.listdir(tmp_dir)
    for tmp_f in dir_list:
        if tmp_f[0:6] == 'temp_p':
            pname = os.path.join(tmp_dir, tmp_f)
            sample_permutation = np.load(pname)
            hname = os.path.join(tmp_dir, tmp_f[0:5] + 'h' + tmp_f[6:len(tmp_f)])
            h_mat = np.load(hname)
            linkage_matrix = keg.update_linkage_matrix(h_mat, sample_permutation, linkage_matrix)

    return linkage_matrix

def save_temporary_cluster(h_matrix, sample_permutation, run_parameters, sequence_number):
    """ save one h_matrix and one permutation in temorary files with sequence_number appended names.

    Args:
        h_matrix: k x permutation size matrix.
        sample_permutation: indices of h_matrix columns permutation.
        run_parameters: parmaeters including the "tmp_directory" name.
        sequence_number: temporary file name suffix.
    """
    tmp_dir = run_parameters["tmp_directory"]
    time_stamp = timestamp_filename('_N', str(sequence_number), run_parameters)
    hname = os.path.join(tmp_dir, 'temp_h'+time_stamp)
    cluster_id = np.argmax(h_matrix, 0)
    cluster_id.dump(hname)
    pname = os.path.join(tmp_dir, 'temp_p'+time_stamp)
    sample_permutation.dump(pname)

    return
    
def form_consensus_matrix_graphic(consensus_matrix, k=3):
    ''' use K-means to reorder the consensus matrix for graphic display.

    Args:
        consensus_matrix: calculated consensus matrix in samples x samples order.
        k: number of clusters estimate (inner diminsion k of factored h_matrix).

    Returns:
        cc_cm: consensus_matrix with rows and columns in K-means sort order.
    '''
    cc_cm = consensus_matrix.copy()
    labels = keg.perform_kmeans(consensus_matrix, k)
    sorted_labels = np.argsort(labels)
    cc_cm = cc_cm[sorted_labels[:, None], sorted_labels]

    return cc_cm
    
def display_clusters(consensus_matrix):
    ''' graphic display the consensus matrix.

    Args:
         consenus matrix: usually a smallish square matrix.
    '''
    methods = [None, 'none', 'nearest', 'bilinear', 'bicubic', 'spline16',
               'spline36', 'hanning', 'hamming', 'hermite', 'kaiser', 'quadric',
               'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos']
    grid = consensus_matrix
    fig, axes = plt.subplots(3, 6, figsize=(12, 6),
                             subplot_kw={'xticks': [], 'yticks': []})
    fig.subplots_adjust(hspace=0.3, wspace=0.05)
    for ax_n, interp_method in zip(axes.flat, methods):
        ax_n.imshow(grid, interpolation=interp_method)
        ax_n.set_title(interp_method)
    plt.show()

    return
    
def save_consensus_cluster_result(consensus_matrix, sample_names, labels, run_parameters):
    """ write the results of network based nmf consensus clustering to output files.

    Args:
        consensus_matrix: sample_names x labels - symmetric consensus matrix.
        sample_names: spreadsheet column names.
        labels: cluster assignments for column names (or consensus matrix).
        run_parameters: dictionary with "results_directory" key.
    """
    write_consensus_matrix(consensus_matrix, sample_names, labels, run_parameters)
    save_clusters(sample_names, labels, run_parameters)

    return

def write_consensus_matrix(consensus_matrix, sample_names, labels, run_parameters):
    """ write the consensus matrix as a dataframe with sample_names column lablels
        and cluster labels as row labels.

    Args:
        consensus_matrix: sample_names x sample_names numerical matrix.
        sample_names: data identifiers for column names.
        labels: cluster numbers for row names.
        run_parameters: path to write to consensus_data file (run_parameters["results_directory"]).
    """
    if int(run_parameters["use_now_name"]) != 0:
        file_name = os.path.join(
            run_parameters["results_directory"], timestamp_filename('consensus_data', 'df'))
    else:
        file_name = os.path.join(run_parameters["results_directory"], 'consensus_data.df')
    out_df = pd.DataFrame(data=consensus_matrix, columns=sample_names, index=labels)
    out_df.to_csv(file_name, sep='\t')

    return

def save_clusters(sample_names, labels, run_parameters):
    """ wtite .tsv file that assings a cluster number label to the sample_names.

    Args:
        sample_names: (unique) data identifiers.
        labels: cluster number assignments.
        run_parameters: write path (run_parameters["results_directory"]).
    """
    if int(run_parameters["use_now_name"]) != 0:
        file_name = os.path.join(
            run_parameters["results_directory"], timestamp_filename('labels_data', 'tsv'))
    else:
        file_name = os.path.join(run_parameters["results_directory"], 'labels_data.tsv')

    df_tmp = pd.DataFrame(data=labels, index=sample_names)
    df_tmp.to_csv(file_name, sep='\t', header=None)

    return
    
def timestamp_filename(name_base, name_extension, run_parameters=None):
    """ insert a time stamp into the filename_ before .extension.

    Args:
        name_base: file name first part - may include directory path.
        name_extension: file extension without a period.
        run_parameters: run_parameters['use_now_name'] (between 1 and 1,000,000)

    Returns:
        time_stamped_file_name: concatenation of time-stamp between inputs.
    """
    dt_max = 1e6
    dt_min = 1
    if run_parameters is None:
        nstr = time.strftime("%a_%d_%b_%Y_%H_%M_%S", time.localtime())
    else:
        time_step = min(max(int(run_parameters['use_now_name']), dt_min), dt_max)
        nstr = np.str_(int(time.time() * time_step))

    time_stamped_file_name = name_base + '_' + nstr + '.' + name_extension

    return time_stamped_file_name
    
def echo_input(network_mat, spreadsheet_mat, run_parameters):
    ''' command line display data: network and spreadsheet matrices and run parameters.

    Args:
         network_mat: gene-gene network matrix.
         spreadsheet_mat: genes x samples user spreadsheet data matrix.
         run_parameters: run parameters dictionary.
    '''
    net_rows = network_mat.shape[0]
    net_cols = network_mat.shape[1]
    usr_rows = spreadsheet_mat.shape[0]
    usr_cols = spreadsheet_mat.shape[1]
    print('\nMethod: {}'.format(run_parameters['method']))
    date_frm = "Local: %a, %d %b %Y %H:%M:%S"
    print('Data Loaded:\t{}'.format(time.strftime(date_frm, time.localtime())))
    print('\nnetwork_file_name: {}'.format(run_parameters['network_file_name']))
    print('network    matrix {} x {}'.format(net_rows, net_cols))
    print('\nsamples_file_name: {}'.format(run_parameters['samples_file_name']))
    print('spread sheet matrix {} x {}\n'.format(usr_rows, usr_cols))
    print('\nAll run parameters as received:\n')
    display_run_parameters(run_parameters)

    return

def display_run_parameters(run_parameters):
    """ command line display the run parameters dictionary.

    Args:
        run_parameters: dictionary of run parameters.
    """
    for fielap_dag_n in run_parameters:
        print('{} : {}'.format(fielap_dag_n, run_parameters[fielap_dag_n]))
    print('\n')

    return

def generate_run_file(run_parameters=None, file_name='run_file'):
    """ write a parameter set dictionary to a text file for editing.

    Args:
        file_name: file name (will be written as plain text).
    """
    if run_parameters is None:
        run_parameters = {
            "method":"cc_net_cluster_nmf",
            "k":4,
            "number_of_bootstraps":5,
            "percent_sample":0.8,
            "restart_probability":0.7,
            "number_of_iteriations_in_rwr":100,
            "it_max":2000,
            "h_clust_eq_limit":200,
            "obj_fcn_chk_freq":50,
            "restart_tolerance":1e-4,
            'lmbda':1400,
            "network_file_name":"network_file_name",
            "samples_file_name":"samples_file_name",
            "tmp_directory":"tmp",
            "results_directory":"results",
            "use_now_name":1,
            "verbose":1,
            "display_clusters":1}

    with open(file_name, 'w') as yaml_file:
        yaml_file.write(yaml.dump(run_parameters, default_flow_style=False))

    return

def compare_cluster_labels(x_file_name, y_file_name):
    """ compare .tsv files to see if the same names are in the same clusters
    """
    x_df = pd.read_csv(x_file_name, header=None, names=None, delimiter='\t', usecols=[0, 1])
    x_df.columns = ['x_ID', 'x_cluster_n']
    x_df = x_df.sort_values('x_ID', ascending=0)
    y_df = pd.read_csv(y_file_name, header=None, names=None, delimiter='\t', usecols=[0, 1])
    y_df.columns = ['y_ID', 'y_cluster_n']
    y_df = y_df.sort_values('y_ID', ascending=0)
    frames = [x_df, y_df]
    comp_df = pd.concat(frames, axis=1)
    x_lbl_arr = np.array(comp_df['x_cluster_n'])
    c_max = x_lbl_arr.max()
    c_min = x_lbl_arr.min()
    
    for cn in range(c_min, c_max + 1):
        tmp = comp_df.loc[(comp_df['x_cluster_n'] == cn)]
        n_tot = tmp.shape[0]
        n = most_common_n(tmp['y_cluster_n'])
        n_not = sum(tmp['y_cluster_n'] != n)
        print('cn={}, n={}\t{} of {} agree, {} disagree'.format(cn, n, n_tot - n_not, n_tot, n_not))

    return comp_df

def most_common_n(n_arr):
    n_arr = np.array(n_arr)
    mx = max(n_arr)
    mn = min(n_arr)
    nc = mn
    for t in range(mn, mx+1):
        if sum(n_arr == t) > sum(n_arr == nc):
            nc = t
    return nc