# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 14:47:45 2016

@author: The Gene Sets Characterization dev team

"""
#import multiprocessing

# Number of processes to be executed in parallel
#number_of_processes = multiprocessing.cpu_count()
#print("Using parallelism {}".format(number_of_processes))

def nmf(run_parameters):
    '''nmf clustering'''
    from sample_clustering_toolbox import run_nmf
    run_nmf(run_parameters) 

def cc_nmf(run_parameters):
    '''kmeans consensus clustering of the nmf-based clusters'''
    from sample_clustering_toolbox import run_cc_nmf
    print('debug tracking: cc_nmf called:')
    print(run_parameters)
    run_cc_nmf(run_parameters)

def net_nmf(run_parameters):
    '''net-nmf clustering "'''
    from sample_clustering_toolbox import run_net_nmf
    run_net_nmf(run_parameters)

def cc_net_nmf(run_parameters):
    '''kmeans consensus clustering of the net-nmf-based clusters'''
    from sample_clustering_toolbox import run_cc_net_nmf
    print('debug tracking: cc_net_nmf called:')
    print(run_parameters)
    run_cc_net_nmf(run_parameters)

SELECT = {
    "cluster_nmf":nmf,
    "cc_cluster_nmf":cc_nmf,
    "net_cluster_nmf":net_nmf,
    "cc_net_cluster_nmf":cc_net_nmf}

def main():
    import sys
    from knpackage.toolbox import get_run_directory_and_file
    from knpackage.toolbox import get_run_parameters
    
    run_directory, run_file = get_run_directory_and_file(sys.argv)
    run_parameters = get_run_parameters(run_directory, run_file)
    SELECT[run_parameters["method"]](run_parameters)

if __name__ == "__main__":
    main()
