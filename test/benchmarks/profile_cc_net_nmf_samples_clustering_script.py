"""
network_file_name: /Users/del/AllCodeBigData/KnowEnG_tbx_test/input_data/final_clean_4col.edge
samples_file_name: /Users/del/AllCodeBigData/KnowEnG_tbx_test/input_data/final_clean_full_matrix.df
results_directory: /Users/del/AllCodeBigData/KnowEnG_tbx_test/results
tmp_directory: /Users/del/AllCodeBigData/KnowEnG_tbx_test/tmp
method: cc_net_cluster_nmf
k: "3"
number_of_iteriations_in_rwr: "100"
obj_fcn_chk_freq: "50"
it_max: "10000"
h_clust_eq_limit: "200"
restart_tolerance: "0.0001"
lmbda: "1400"
percent_sample: "0.8"
number_of_bootstraps: "5"
display_clusters: "1"
restart_probability: "0.7"
verbose: "1"
use_now_name: "1000000"
method1: cluster_nmf
method2: cc_cluster_nmf
method3: net_cluster_nmf
method4: cc_net_cluster_nmf
"""

from knpackage.toolbox import get_run_parameters
from sample_clustering_toolbox import run_cc_net_nmf

run_directory = '/Users/lanier4/BigDataTank/nbs_run'
run_file =  'cc_net_cluster_nmf_run_file.yml'
run_parameters = get_run_parameters(run_directory, run_file)
run_parameters['display_clusters'] = 0
run_parameters['number_of_bootstraps'] = 1

run_cc_net_nmf(run_parameters)