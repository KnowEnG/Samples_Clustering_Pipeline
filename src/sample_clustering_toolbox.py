# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 16:08:25 2016

@author: del
@author: The Gene Sets Characterization dev team

"""
import os
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import knpackage.toolbox as kn


def run_nmf(run_parameters):
    """ wrapper: call sequence to perform non-negative matrix factorization and write results.

    Args:
        run_parameters: parameter set dictionary.
    """
    spreadsheet_df = kn.get_spreadsheet_df(run_parameters['spreadsheet_name_full_path'])
    spreadsheet_mat = spreadsheet_df.as_matrix()
    spreadsheet_mat = kn.get_quantile_norm_matrix(spreadsheet_mat)

    h_mat = kn.perform_nmf(spreadsheet_mat, run_parameters)

    linkage_matrix = np.zeros((spreadsheet_mat.shape[1], spreadsheet_mat.shape[1]))
    sample_perm = np.arange(0, spreadsheet_mat.shape[1])
    linkage_matrix = kn.update_linkage_matrix(h_mat, sample_perm, linkage_matrix)
    labels = kn.perform_kmeans(linkage_matrix, int(run_parameters['number_of_clusters']))

    sample_names = spreadsheet_df.columns
    save_final_samples_clustering(sample_names, labels, run_parameters)

    if int(run_parameters['display_clusters']) != 0:
        con_mat_image = form_consensus_matrix_graphic(linkage_matrix, int(run_parameters['number_of_clusters']))
        display_clusters(con_mat_image)

    return


def run_cc_nmf(run_parameters):
    """ wrapper: call sequence to perform non-negative matrix factorization with
        consensus clustering and write results.

    Args:
        run_parameters: parameter set dictionary.
    """
    tmp_dir = 'tmp_cc_nmf'
    run_parameters["tmp_directory"] = kn.create_dir(
        run_parameters["run_directory"], tmp_dir)

    spreadsheet_df = kn.get_spreadsheet_df(run_parameters['spreadsheet_name_full_path'])
    spreadsheet_mat = spreadsheet_df.as_matrix()
    spreadsheet_mat = kn.get_quantile_norm_matrix(spreadsheet_mat)

    find_and_save_nmf_clusters(spreadsheet_mat, run_parameters)

    linkage_matrix = np.zeros((spreadsheet_mat.shape[1], spreadsheet_mat.shape[1]))
    indicator_matrix = linkage_matrix.copy()
    consensus_matrix = form_consensus_matrix(run_parameters, linkage_matrix, indicator_matrix)
    labels = kn.perform_kmeans(consensus_matrix, int(run_parameters['number_of_clusters']))

    sample_names = spreadsheet_df.columns
    save_consensus_clustering(consensus_matrix, sample_names, labels, run_parameters)
    save_final_samples_clustering(sample_names, labels, run_parameters)

    kn.remove_dir(run_parameters["tmp_directory"])

    if int(run_parameters['display_clusters']) != 0:
        display_clusters(form_consensus_matrix_graphic(consensus_matrix, int(run_parameters['number_of_clusters'])))

    return


def run_net_nmf(run_parameters):
    """ wrapper: call sequence to perform network based stratification and write results.

    Args:
        run_parameters: parameter set dictionary.
    """
    spreadsheet_df = kn.get_spreadsheet_df(run_parameters['spreadsheet_name_full_path'])
    network_df = kn.get_network_df(run_parameters['gg_network_name_full_path'])

    node_1_names, node_2_names = kn.extract_network_node_names(network_df)
    unique_gene_names = kn.find_unique_node_names(node_1_names, node_2_names)
    genes_lookup_table = kn.create_node_names_dict(unique_gene_names)

    network_df = kn.map_node_names_to_index(network_df, genes_lookup_table, 'node_1')
    network_df = kn.map_node_names_to_index(network_df, genes_lookup_table, 'node_2')

    network_df = kn.symmetrize_df(network_df)
    network_mat = kn.convert_network_df_to_sparse(
        network_df, len(unique_gene_names), len(unique_gene_names))

    network_mat = kn.normalize_sparse_mat_by_diagonal(network_mat)
    lap_diag, lap_pos = kn.form_network_laplacian_matrix(network_mat)

    spreadsheet_df = kn.update_spreadsheet_df(spreadsheet_df, unique_gene_names)
    spreadsheet_mat = spreadsheet_df.as_matrix()
    sample_names = spreadsheet_df.columns

    sample_smooth, iterations = kn.smooth_matrix_with_rwr(
        spreadsheet_mat, network_mat, run_parameters)
    sample_quantile_norm = kn.get_quantile_norm_matrix(sample_smooth)
    h_mat = kn.perform_net_nmf(sample_quantile_norm, lap_pos, lap_diag, run_parameters)

    linkage_matrix = np.zeros((spreadsheet_mat.shape[1], spreadsheet_mat.shape[1]))
    sample_perm = np.arange(0, spreadsheet_mat.shape[1])
    linkage_matrix = kn.update_linkage_matrix(h_mat, sample_perm, linkage_matrix)
    labels = kn.perform_kmeans(linkage_matrix, int(run_parameters["number_of_clusters"]))

    save_final_samples_clustering(sample_names, labels, run_parameters)

    if int(run_parameters['display_clusters']) != 0:
        display_clusters(form_consensus_matrix_graphic(linkage_matrix, int(run_parameters['number_of_clusters'])))

    return


def run_cc_net_nmf(run_parameters):
    """ wrapper: call sequence to perform network based stratification with consensus clustering
        and write results.

    Args:
        run_parameters: parameter set dictionary.
    """
    tmp_dir = 'tmp_cc_net_nmf'
    run_parameters["tmp_directory"] = kn.create_dir(
        run_parameters["run_directory"], tmp_dir)

    spreadsheet_df = kn.get_spreadsheet_df(run_parameters['spreadsheet_name_full_path'])
    network_df = kn.get_network_df(run_parameters['gg_network_name_full_path'])

    node_1_names, node_2_names = kn.extract_network_node_names(network_df)
    unique_gene_names = kn.find_unique_node_names(node_1_names, node_2_names)
    genes_lookup_table = kn.create_node_names_dict(unique_gene_names)

    network_df = kn.map_node_names_to_index(network_df, genes_lookup_table, 'node_1')
    network_df = kn.map_node_names_to_index(network_df, genes_lookup_table, 'node_2')

    network_df = kn.symmetrize_df(network_df)
    # network_mat = convert_df_to_sparse(network_df, len(unique_gene_names))
    network_mat = kn.convert_network_df_to_sparse(
        network_df, len(unique_gene_names), len(unique_gene_names))

    network_mat = kn.normalize_sparse_mat_by_diagonal(network_mat)
    lap_diag, lap_pos = kn.form_network_laplacian_matrix(network_mat)

    spreadsheet_df = kn.update_spreadsheet_df(spreadsheet_df, unique_gene_names)
    spreadsheet_mat = spreadsheet_df.as_matrix()
    sample_names = spreadsheet_df.columns

    find_and_save_net_nmf_clusters(network_mat, spreadsheet_mat, lap_diag, lap_pos, run_parameters)

    linkage_matrix = np.zeros((spreadsheet_mat.shape[1], spreadsheet_mat.shape[1]))
    indicator_matrix = linkage_matrix.copy()
    consensus_matrix = form_consensus_matrix(run_parameters, linkage_matrix, indicator_matrix)
    labels = kn.perform_kmeans(consensus_matrix, int(run_parameters['number_of_clusters']))

    save_consensus_clustering(consensus_matrix, sample_names, labels, run_parameters)
    save_final_samples_clustering(sample_names, labels, run_parameters)

    kn.remove_dir(run_parameters["tmp_directory"])

    if int(run_parameters['display_clusters']) != 0:
        display_clusters(form_consensus_matrix_graphic(consensus_matrix, int(run_parameters['number_of_clusters'])))

    return


def find_and_save_net_nmf_clusters(network_mat, spreadsheet_mat, lap_dag, lap_val, run_parameters):
    """ central loop: compute components for the consensus matrix from the input
        network and spreadsheet matrices and save them to temp files.

    Args:
        network_mat: genes x genes symmetric matrix.
        spreadsheet_mat: genes x samples matrix.
        lap_dag: laplacian matrix component, L = lap_dag - lap_val.
        lap_val: laplacian matrix component, L = lap_dag - lap_val.
        run_parameters: dictionay of run-time parameters.
    """
    for sample in range(0, int(run_parameters["number_of_bootstraps"])):
        sample_random, sample_permutation = kn.sample_a_matrix(
            spreadsheet_mat, np.float64(run_parameters["rows_sampling_fraction"]),
            np.float64(run_parameters["cols_sampling_fraction"]))
        sample_smooth, iterations = \
            kn.smooth_matrix_with_rwr(sample_random, network_mat, run_parameters)

        print("{} of {}: iterations = {}".format(sample + 1, run_parameters["number_of_bootstraps"], iterations))

        sample_quantile_norm = kn.get_quantile_norm_matrix(sample_smooth)
        h_mat = kn.perform_net_nmf(sample_quantile_norm, lap_val, lap_dag, run_parameters)

        save_a_clustering_to_tmp(h_mat, sample_permutation, run_parameters, sample)

    return


def find_and_save_nmf_clusters(spreadsheet_mat, run_parameters):
    """ central loop: compute components for the consensus matrix by
        non-negative matrix factorization.

    Args:
        spreadsheet_mat: genes x samples matrix.
        run_parameters: dictionay of run-time parameters.
    """
    for sample in range(0, int(run_parameters["number_of_bootstraps"])):
        sample_random, sample_permutation = kn.sample_a_matrix(
            spreadsheet_mat, np.float64(run_parameters["rows_sampling_fraction"]),
            np.float64(run_parameters["cols_sampling_fraction"]))

        h_mat = kn.perform_nmf(sample_random, run_parameters)
        save_a_clustering_to_tmp(h_mat, sample_permutation, run_parameters, sample)

        print('nmf {} of {}'.format(sample + 1, run_parameters["number_of_bootstraps"]))

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
    indicator_matrix = get_indicator_matrix(run_parameters, indicator_matrix)
    linkage_matrix = get_linkage_matrix(run_parameters, linkage_matrix)
    consensus_matrix = linkage_matrix / np.maximum(indicator_matrix, 1)

    return consensus_matrix


def get_indicator_matrix(run_parameters, indicator_matrix):
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
            indicator_matrix = kn.update_indicator_matrix(sample_permutation, indicator_matrix)

    return indicator_matrix


def get_linkage_matrix(run_parameters, linkage_matrix):
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
            linkage_matrix = kn.update_linkage_matrix(h_mat, sample_permutation, linkage_matrix)

    return linkage_matrix


def save_a_clustering_to_tmp(h_matrix, sample_permutation, run_parameters, sequence_number):
    """ save one h_matrix and one permutation in temorary files with sequence_number appended names.

    Args:
        h_matrix: k x permutation size matrix.
        sample_permutation: indices of h_matrix columns permutation.
        run_parameters: parmaeters including the "tmp_directory" name.
        sequence_number: temporary file name suffix.
    """
    tmp_dir = run_parameters["tmp_directory"]
    #time_stamp = timestamp_filename('_N', str(sequence_number), run_parameters)
    time_stamp = kn.create_timestamped_filename('_N' + str(sequence_number), name_extension=None, precision=1e12)

    hname = os.path.join(tmp_dir, 'temp_h' + time_stamp)
    pname = os.path.join(tmp_dir, 'temp_p' + time_stamp)

    cluster_id = np.argmax(h_matrix, 0)
    with open(hname, 'wb') as fh0:
        cluster_id.dump(fh0)
    with open(pname, 'wb') as fh1:
        sample_permutation.dump(fh1)

    return


def form_consensus_matrix_graphic(consensus_matrix, k=3):
    """ use K-means to reorder the consensus matrix for graphic display.

    Args:
        consensus_matrix: calculated consensus matrix in samples x samples order.
        k: number of clusters estimate (inner diminsion k of factored h_matrix).

    Returns:
        cc_cm: consensus_matrix with rows and columns in K-means sort order.
    """
    cc_cm = consensus_matrix.copy()
    labels = kn.perform_kmeans(consensus_matrix, k)
    sorted_labels = np.argsort(labels)
    cc_cm = cc_cm[sorted_labels[:, None], sorted_labels]

    return cc_cm


def display_clusters(consensus_matrix):
    """ graphic display the consensus matrix.

    Args:
         consenus matrix: usually a smallish square matrix.
    """
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


def save_consensus_clustering(consensus_matrix, sample_names, labels, run_parameters):
    """ write the consensus matrix as a dataframe with sample_names column lablels
        and cluster labels as row labels.

    Args:
        consensus_matrix: sample_names x sample_names numerical matrix.
        sample_names: data identifiers for column names.
        labels: cluster numbers for row names.
        run_parameters: path to write to consensus_data file (run_parameters["results_directory"]).
    """
    file_name = os.path.join(run_parameters["results_directory"], kn.create_timestamped_filename('consensus_data', 'df'))
    out_df = pd.DataFrame(data=consensus_matrix, columns=sample_names, index=labels)
    out_df.to_csv(file_name, sep='\t')

    return


def save_final_samples_clustering(sample_names, labels, run_parameters):
    """ wtite .tsv file that assings a cluster number label to the sample_names.

    Args:
        sample_names: (unique) data identifiers.
        labels: cluster number assignments.
        run_parameters: write path (run_parameters["results_directory"]).
    """
    file_name = os.path.join(run_parameters["results_directory"], kn.create_timestamped_filename('labels_data', 'tsv'))
    df_tmp = kn.create_df_with_sample_labels(sample_names, labels)
    df_tmp.to_csv(file_name, sep='\t', header=None)

    return