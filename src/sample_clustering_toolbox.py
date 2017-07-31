"""
@author: The KnowEnG dev team
"""
import os
import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score
import knpackage.toolbox as kn
import knpackage.distributed_computing_utils as dstutil

import clustering_eval_toolbox as cluster_eval

def run_nmf(run_parameters):
    """ wrapper: call sequence to perform non-negative matrix factorization and write results.

    Args:
        run_parameters: parameter set dictionary.
    """
    number_of_clusters = run_parameters['number_of_clusters']
    spreadsheet_name_full_path = run_parameters['spreadsheet_name_full_path']

    spreadsheet_df = kn.get_spreadsheet_df(spreadsheet_name_full_path)
    spreadsheet_mat = spreadsheet_df.as_matrix()
    spreadsheet_mat = kn.get_quantile_norm_matrix(spreadsheet_mat)

    h_mat = kn.perform_nmf(spreadsheet_mat, run_parameters)

    linkage_matrix = np.zeros((spreadsheet_mat.shape[1], spreadsheet_mat.shape[1]))
    sample_perm = np.arange(0, spreadsheet_mat.shape[1])
    linkage_matrix = kn.update_linkage_matrix(h_mat, sample_perm, linkage_matrix)
    labels = kn.perform_kmeans(linkage_matrix, number_of_clusters)

    sample_names = spreadsheet_df.columns
    save_consensus_clustering(linkage_matrix, sample_names, labels, run_parameters)
    save_final_samples_clustering(sample_names, labels, run_parameters)
    save_spreadsheet_and_variance_heatmap(spreadsheet_df, labels, run_parameters)


def run_net_nmf(run_parameters):
    """ wrapper: call sequence to perform network based stratification and write results.

    Args:
        run_parameters: parameter set dictionary.
    """
    number_of_clusters = run_parameters['number_of_clusters']
    gg_network_name_full_path = run_parameters['gg_network_name_full_path']
    spreadsheet_name_full_path = run_parameters['spreadsheet_name_full_path']

    network_mat, unique_gene_names = kn.get_sparse_network_matrix(gg_network_name_full_path)
    network_mat = kn.normalize_sparse_mat_by_diagonal(network_mat)
    lap_diag, lap_pos = kn.form_network_laplacian_matrix(network_mat)

    spreadsheet_df = kn.get_spreadsheet_df(spreadsheet_name_full_path)
    spreadsheet_df = kn.update_spreadsheet_df(spreadsheet_df, unique_gene_names)
    sample_names = spreadsheet_df.columns
    spreadsheet_mat = spreadsheet_df.as_matrix()

    spreadsheet_mat, iterations = kn.smooth_matrix_with_rwr(spreadsheet_mat, network_mat, run_parameters)
    spreadsheet_mat = kn.get_quantile_norm_matrix(spreadsheet_mat)
    h_mat = kn.perform_net_nmf(spreadsheet_mat, lap_pos, lap_diag, run_parameters)

    linkage_matrix = np.zeros((spreadsheet_mat.shape[1], spreadsheet_mat.shape[1]))
    sample_perm = np.arange(0, spreadsheet_mat.shape[1])
    linkage_matrix = kn.update_linkage_matrix(h_mat, sample_perm, linkage_matrix)
    labels = kn.perform_kmeans(linkage_matrix, number_of_clusters)

    save_consensus_clustering(linkage_matrix, sample_names, labels, run_parameters)
    save_final_samples_clustering(sample_names, labels, run_parameters)
    save_spreadsheet_and_variance_heatmap(spreadsheet_df, labels, run_parameters, network_mat)


def run_cc_nmf(run_parameters):
    """ wrapper: call sequence to perform non-negative matrix factorization with
        consensus clustering and write results.

    Args:
        run_parameters: parameter set dictionary.
    """
    tmp_dir = 'tmp_cc_nmf'
    run_parameters = update_tmp_directory(run_parameters, tmp_dir)

    processing_method = run_parameters['processing_method']
    number_of_bootstraps = run_parameters['number_of_bootstraps']
    number_of_clusters = run_parameters['number_of_clusters']
    spreadsheet_name_full_path = run_parameters['spreadsheet_name_full_path']

    spreadsheet_df = kn.get_spreadsheet_df(spreadsheet_name_full_path)
    spreadsheet_mat = spreadsheet_df.as_matrix()
    spreadsheet_mat = kn.get_quantile_norm_matrix(spreadsheet_mat)
    number_of_samples = spreadsheet_mat.shape[1]

    if processing_method == 'serial':
        for sample in range(0, number_of_bootstraps):
            run_cc_nmf_clusters_worker(spreadsheet_mat, run_parameters, sample)

    elif processing_method == 'parallel':
        find_and_save_cc_nmf_clusters_parallel(spreadsheet_mat, run_parameters, number_of_bootstraps)

    elif processing_method == 'distribute':
        func_args = [spreadsheet_mat, run_parameters]
        dependency_list = [run_cc_nmf_clusters_worker, save_a_clustering_to_tmp, dstutil.determine_parallelism_locally]
        dstutil.execute_distribute_computing_job(run_parameters['cluster_ip_address'],
                                                 number_of_bootstraps,
                                                 func_args,
                                                 find_and_save_cc_nmf_clusters_parallel,
                                                 dependency_list)
    else:
        raise ValueError('processing_method contains bad value.')

    consensus_matrix = form_consensus_matrix(run_parameters, number_of_samples)
    labels = kn.perform_kmeans(consensus_matrix, number_of_clusters)

    sample_names = spreadsheet_df.columns
    save_consensus_clustering(consensus_matrix, sample_names, labels, run_parameters)
    save_final_samples_clustering(sample_names, labels, run_parameters)
    save_spreadsheet_and_variance_heatmap(spreadsheet_df, labels, run_parameters)

    kn.remove_dir(run_parameters["tmp_directory"])


def run_cc_net_nmf(run_parameters):
    """ wrapper: call sequence to perform network based stratification with consensus clustering
        and write results.

    Args:
        run_parameters: parameter set dictionary.
    """
    tmp_dir = 'tmp_cc_net_nmf'
    run_parameters = update_tmp_directory(run_parameters, tmp_dir)

    processing_method = run_parameters['processing_method']
    number_of_bootstraps = run_parameters['number_of_bootstraps']
    number_of_clusters = run_parameters['number_of_clusters']
    gg_network_name_full_path = run_parameters['gg_network_name_full_path']
    spreadsheet_name_full_path = run_parameters['spreadsheet_name_full_path']

    network_mat, unique_gene_names = kn.get_sparse_network_matrix(gg_network_name_full_path)
    network_mat = kn.normalize_sparse_mat_by_diagonal(network_mat)
    lap_diag, lap_pos = kn.form_network_laplacian_matrix(network_mat)

    spreadsheet_df = kn.get_spreadsheet_df(spreadsheet_name_full_path)
    spreadsheet_df = kn.update_spreadsheet_df(spreadsheet_df, unique_gene_names)
    spreadsheet_mat = spreadsheet_df.as_matrix()
    number_of_samples = spreadsheet_mat.shape[1]
    sample_names = spreadsheet_df.columns

    if processing_method == 'serial':
        for sample in range(0, number_of_bootstraps):
            run_cc_net_nmf_clusters_worker(network_mat, spreadsheet_mat, lap_diag, lap_pos, run_parameters, sample)

    elif processing_method == 'parallel':
        find_and_save_cc_net_nmf_clusters_parallel(network_mat, spreadsheet_mat, lap_diag, lap_pos, run_parameters, number_of_bootstraps)

    elif processing_method == 'distribute':
        func_args = [network_mat, spreadsheet_mat, lap_diag, lap_pos, run_parameters]
        dependency_list = [run_cc_net_nmf_clusters_worker, save_a_clustering_to_tmp, dstutil.determine_parallelism_locally]
        dstutil.execute_distribute_computing_job(run_parameters['cluster_ip_address'],
                                                 number_of_bootstraps,
                                                 func_args,
                                                 find_and_save_cc_net_nmf_clusters_parallel,
                                                 dependency_list)
    else:
        raise ValueError('processing_method contains bad value.')

    consensus_matrix = form_consensus_matrix(run_parameters, number_of_samples)
    labels = kn.perform_kmeans(consensus_matrix, number_of_clusters)

    save_consensus_clustering(consensus_matrix, sample_names, labels, run_parameters)
    save_final_samples_clustering(sample_names, labels, run_parameters)
    save_spreadsheet_and_variance_heatmap(spreadsheet_df, labels, run_parameters, network_mat)

    kn.remove_dir(run_parameters["tmp_directory"])


def find_and_save_cc_nmf_clusters_parallel(spreadsheet_mat, run_parameters, local_parallelism):
    """ central loop: compute components for the consensus matrix by
        non-negative matrix factorization.

    Args:
        spreadsheet_mat: genes x samples matrix.
        run_parameters: dictionary of run-time parameters.
        number_of_cpus: number of processes to be running in parallel
    """
    import knpackage.distributed_computing_utils as dstutil

    jobs_id = range(0, local_parallelism)
    zipped_arguments = dstutil.zip_parameters(spreadsheet_mat, run_parameters, jobs_id)
    if 'parallelism' in run_parameters:
        parallelism = dstutil.determine_parallelism_locally(local_parallelism, run_parameters['parallelism'])
    else:
        parallelism = dstutil.determine_parallelism_locally(local_parallelism)
    dstutil.parallelize_processes_locally(run_cc_nmf_clusters_worker, zipped_arguments, parallelism)


def find_and_save_cc_net_nmf_clusters_parallel(network_mat, spreadsheet_mat, lap_diag, lap_pos, run_parameters, local_parallelism):
    """ central loop: compute components for the consensus matrix from the input
        network and spreadsheet matrices and save them to temp files.

    Args:
        network_mat: genes x genes symmetric matrix.
        spreadsheet_mat: genes x samples matrix.
        lap_dag: laplacian matrix component, L = lap_dag - lap_val.
        lap_val: laplacian matrix component, L = lap_dag - lap_val.
        run_parameters: dictionary of run-time parameters.
        number_of_cpus: number of processes to be running in parallel
    """
    import knpackage.distributed_computing_utils as dstutil

    jobs_id = range(0, local_parallelism)
    zipped_arguments = dstutil.zip_parameters(network_mat, spreadsheet_mat, lap_diag, lap_pos, run_parameters, jobs_id)
    if 'parallelism' in run_parameters:
        parallelism = dstutil.determine_parallelism_locally(local_parallelism, run_parameters['parallelism'])
    else:
        parallelism = dstutil.determine_parallelism_locally(local_parallelism)
    dstutil.parallelize_processes_locally(run_cc_net_nmf_clusters_worker, zipped_arguments, parallelism)


def run_cc_nmf_clusters_worker(spreadsheet_mat, run_parameters, sample):
    """Worker to execute nmf_clusters in a single process

    Args:
        spreadsheet_mat: genes x samples matrix.
        run_parameters: dictionary of run-time parameters.
        sample: each loops.

    Returns:
        None

    """
    import knpackage.toolbox as kn
    import numpy as np

    np.random.seed(sample)
    rows_sampling_fraction = run_parameters["rows_sampling_fraction"]
    cols_sampling_fraction = run_parameters["cols_sampling_fraction"]
    spreadsheet_mat, sample_permutation = kn.sample_a_matrix(spreadsheet_mat,
                                                             rows_sampling_fraction, cols_sampling_fraction)
    h_mat = kn.perform_nmf(spreadsheet_mat, run_parameters)
    save_a_clustering_to_tmp(h_mat, sample_permutation, run_parameters, sample)


def run_cc_net_nmf_clusters_worker(network_mat, spreadsheet_mat, lap_dag, lap_val, run_parameters, sample):
    """Worker to execute net_nmf_clusters in a single process

    Args:
        network_mat: genes x genes symmetric matrix.
        spreadsheet_mat: genes x samples matrix.
        lap_dag: laplacian matrix component, L = lap_dag - lap_val.
        lap_val: laplacian matrix component, L = lap_dag - lap_val.
        run_parameters: dictionay of run-time parameters.
        sample: each single loop.

    Returns:
        None
    """
    import knpackage.toolbox as kn
    import numpy as np

    np.random.seed(sample)
    rows_sampling_fraction = run_parameters["rows_sampling_fraction"]
    cols_sampling_fraction = run_parameters["cols_sampling_fraction"]
    spreadsheet_mat, sample_permutation = kn.sample_a_matrix(spreadsheet_mat, rows_sampling_fraction, cols_sampling_fraction)
    
    spreadsheet_mat, iterations = kn.smooth_matrix_with_rwr(spreadsheet_mat, network_mat, run_parameters)
    spreadsheet_mat = kn.get_quantile_norm_matrix(spreadsheet_mat)
    h_mat = kn.perform_net_nmf(spreadsheet_mat, lap_val, lap_dag, run_parameters)

    save_a_clustering_to_tmp(h_mat, sample_permutation, run_parameters, sample)


def save_a_clustering_to_tmp(h_matrix, sample_permutation, run_parameters, sequence_number):
    """ save one h_matrix and one permutation in temorary files with sequence_number appended names.

    Args:
        h_matrix: k x permutation size matrix.
        sample_permutation: indices of h_matrix columns permutation.
        run_parameters: parmaeters including the "tmp_directory" name.
        sequence_number: temporary file name suffix.
    """
    import os
    import numpy as np

    tmp_dir = run_parameters["tmp_directory"]

    os.makedirs(tmp_dir, mode=0o755, exist_ok=True)

    hname = os.path.join(tmp_dir, 'tmp_h_%d'%(sequence_number))
    pname = os.path.join(tmp_dir, 'tmp_p_%d'%(sequence_number))

    cluster_id = np.argmax(h_matrix, 0)
    with open(hname, 'wb') as fh0:
        cluster_id.dump(fh0)
    with open(pname, 'wb') as fh1:
        sample_permutation.dump(fh1)


def form_consensus_matrix(run_parameters, number_of_samples):
    """ compute the consensus matrix from the indicator and linkage matrix inputs
        formed by the bootstrap "temp_*" files.

    Args:
        run_parameters: parameter set dictionary with "tmp_directory" key.
        linkage_matrix: linkage matrix from initialization or previous call.
        indicator_matrix: indicator matrix from initialization or previous call.

    Returns:
        consensus_matrix: (sum of linkage matrices) / (sum of indicator matrices).
    """
    linkage_matrix = np.zeros((number_of_samples, number_of_samples))
    indicator_matrix = linkage_matrix.copy()

    linkage_matrix, indicator_matrix = get_linkage_matrix(run_parameters, linkage_matrix, indicator_matrix)
    consensus_matrix = linkage_matrix / np.maximum(indicator_matrix, 1)

    return consensus_matrix


def get_linkage_matrix(run_parameters, linkage_matrix, indicator_matrix):
    """ read bootstrap temp_h* and temp_p* files, compute and add the linkage_matrix.

    Args:
        run_parameters: parameter set dictionary.
        linkage_matrix: connectivity matrix from initialization or previous call.

    Returns:
        linkage_matrix: summed with "temp_h*" files in run_parameters["tmp_directory"].
    """
    if run_parameters['processing_method'] == 'distribute':
        tmp_dir = os.path.join(run_parameters['cluster_shared_volumn'],
                               os.path.basename(os.path.normpath(run_parameters['tmp_directory'])))
    else:
        tmp_dir = run_parameters["tmp_directory"]
        
    dir_list = os.listdir(tmp_dir)
    for tmp_f in dir_list:
        if tmp_f[0:6] == 'tmp_p_':
            pname = os.path.join(tmp_dir, tmp_f)
            hname = os.path.join(tmp_dir, 'tmp_h_' + tmp_f[6:len(tmp_f)])

            sample_permutation = np.load(pname)
            h_mat = np.load(hname)

            linkage_matrix = kn.update_linkage_matrix(h_mat, sample_permutation, linkage_matrix)
            indicator_matrix = kn.update_indicator_matrix(sample_permutation, indicator_matrix)

    return linkage_matrix, indicator_matrix


def save_spreadsheet_and_variance_heatmap(spreadsheet_df, labels, run_parameters, network_mat=None):
    """ save the full genes by samples spreadsheet as processed or smoothed if network is provided.
        Also save variance in separate file.
    Args:
        spreadsheet_df: the dataframe as processed
        run_parameters: with keys for "results_directory", "method", (optional - "top_number_of_genes")
        network_mat:    (if appropriate) normalized network adjacency matrix used in processing

    Output:
        genes_by_samples_heatmp_{method}_{timestamp}_viz.tsv
        genes_averages_by_cluster_{method}_{timestamp}_viz.tsv
        top_genes_by_cluster_{method}_{timestamp}_download.tsv
    """
    if network_mat is not None:
        sample_smooth, nun = kn.smooth_matrix_with_rwr(spreadsheet_df.as_matrix(), network_mat, run_parameters)
        clusters_df = pd.DataFrame(sample_smooth, index=spreadsheet_df.index.values, columns=spreadsheet_df.columns.values)
    else:
        clusters_df = spreadsheet_df

    clusters_df.to_csv(get_output_file_name(run_parameters, 'genes_by_samples_heatmap', 'viz'), sep='\t')

    cluster_ave_df = pd.DataFrame({i: spreadsheet_df.iloc[:, labels == i].mean(axis=1) for i in np.unique(labels)})
    col_labels = []
    for cluster_number in np.unique(labels):
        col_labels.append('Cluster_%d'%(cluster_number))
    cluster_ave_df.columns = col_labels
    cluster_ave_df.to_csv(get_output_file_name(run_parameters, 'genes_averages_by_cluster', 'viz'), sep='\t')

    clusters_variance_df = pd.DataFrame(clusters_df.var(axis=1), columns=['variance'])
    clusters_variance_df.to_csv(get_output_file_name(run_parameters, 'genes_variance', 'viz'), sep='\t')

    top_number_of_genes_df = pd.DataFrame(data=np.zeros((cluster_ave_df.shape)), columns=cluster_ave_df.columns,
                                          index=cluster_ave_df.index.values)

    top_number_of_genes = run_parameters['top_number_of_genes']
    for sample in top_number_of_genes_df.columns.values:
        top_index = np.argsort(cluster_ave_df[sample].values)[::-1]
        top_number_of_genes_df[sample].iloc[top_index[0:top_number_of_genes]] = 1
    top_number_of_genes_df.to_csv(get_output_file_name(run_parameters, 'top_genes_by_cluster', 'download'), sep='\t')


def save_consensus_clustering(consensus_matrix, sample_names, labels, run_parameters):
    """ write the consensus matrix as a dataframe with sample_names column lablels
        and cluster labels as row labels.

    Args:
        consensus_matrix: sample_names x sample_names numerical matrix.
        sample_names: data identifiers for column names.
        labels: cluster numbers for row names.
        run_parameters: path to write to consensus_data file (run_parameters["results_directory"]).

    Output:
        consensus_matrix_{method}_{timestamp}_viz.tsv
        silhouette_average_{method}_{timestamp}_viz.tsv
    """
    out_df = pd.DataFrame(data=consensus_matrix, columns=sample_names, index=sample_names)
    out_df.to_csv(get_output_file_name(run_parameters, 'consensus_matrix', 'viz'), sep='\t')

    if number_of_labels < 2 or number_of_labels > number_of_samples-1:
       silhouette_average = 1.0
    else:
       silhouette_average = silhouette_score(consensus_matrix, labels)

    silhouette_average = silhouette_score(consensus_matrix, labels)
    silhouette_score_string = 'silhouette number of clusters = %d, corresponding silhouette score = %g' % (
        run_parameters['number_of_clusters'], silhouette_average)

    with open(get_output_file_name(run_parameters, 'silhouette_average', 'viz'), 'w') as fh:
        fh.write(silhouette_score_string)


def save_final_samples_clustering(sample_names, labels, run_parameters):
    """ wtite .tsv file that assings a cluster number label to the sample_names.

    Args:
        sample_names: (unique) data identifiers.
        labels: cluster number assignments.
        run_parameters: write path (run_parameters["results_directory"]).

    Output:
        samples_labeled_by_cluster_{method}_{timestamp}_viz.tsv
        phenotypes_labeled_by_cluster_{method}_{timestamp}_viz.tsv
    """
    cluster_labels_df = kn.create_df_with_sample_labels(sample_names, labels)
    cluster_mapping_full_path = get_output_file_name(run_parameters, 'samples_label_by_cluster', 'viz')
    cluster_labels_df.to_csv(cluster_mapping_full_path, sep='\t', header=None)

    if 'phenotype_name_full_path' in run_parameters.keys():
        run_parameters['cluster_mapping_full_path'] = cluster_mapping_full_path
        cluster_eval.clustering_evaluation(run_parameters)

def get_output_file_name(run_parameters, prefix_string, suffix_string='', type_suffix='tsv'):
    """ get the full directory / filename for writing
    Args:
        run_parameters: dictionary with keys: "results_directory", "method" and "correlation_measure"
        prefix_string:  the first letters of the ouput file name
        suffix_string:  the last letters of the output file name before '.tsv'

    Returns:
        output_file_name:   full file and directory name suitable for file writing
    """
    output_file_name = os.path.join(run_parameters["results_directory"], prefix_string + '_' + run_parameters['method'])
    output_file_name = kn.create_timestamped_filename(output_file_name) + '_' + suffix_string + '.' + type_suffix

    return output_file_name


def update_tmp_directory(run_parameters, tmp_dir):
    ''' Update tmp_directory value in rum_parameters dictionary

    Args:
        run_parameters: run_parameters as the dictionary config
        tmp_dir: temporary directory prefix subjected to different functions

    Returns:
        run_parameters: an updated run_parameters

    '''
    if (run_parameters['processing_method'] == 'distribute'):
        run_parameters["tmp_directory"] = kn.create_dir(run_parameters['cluster_shared_volumn'], tmp_dir)
    else:
        run_parameters["tmp_directory"] = kn.create_dir(run_parameters["run_directory"], tmp_dir)

    return run_parameters


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
