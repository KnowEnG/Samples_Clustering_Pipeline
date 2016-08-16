# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 14:47:45 2016

@author: The Gene Sets Characterization dev team

"""
import time

def disp_run_parameters(run_pars):
    for p in run_pars:
        print('{} : {}'.format(p, run_pars[p]))

def nmf(run_parameters):
    '''nmf clustering'''
    from sample_clustering_toolbox import run_nmf
    t0 = time.time()
    run_nmf(run_parameters)
    print('run_nmf total time: {}'.format(time.time() - t0))
    print('\nWith parameter set:\n')
    disp_run_parameters(run_parameters)

def cc_nmf(run_parameters):
    '''kmeans consensus clustering of the nmf-based clusters'''
    from sample_clustering_toolbox import run_cc_nmf
    t0 = time.time()
    run_cc_nmf(run_parameters)
    print('run_nmf total time: {}'.format(time.time() - t0))
    print('\nWith parameter set:\n')
    disp_run_parameters(run_parameters)

def net_nmf(run_parameters):
    '''net-nmf clustering "'''
    from sample_clustering_toolbox import run_net_nmf
    t0 = time.time()
    run_net_nmf(run_parameters)
    print('run_nmf total time: {}'.format(time.time() - t0))
    print('\nWith parameter set:\n')
    disp_run_parameters(run_parameters)

def cc_net_nmf(run_parameters):
    '''kmeans consensus clustering of the net-nmf-based clusters'''
    from sample_clustering_toolbox import run_cc_net_nmf
    t0 = time.time()
    run_cc_net_nmf(run_parameters)
    print('run_nmf total time: {}'.format(time.time() - t0))
    print('\nWith parameter set:\n')
    disp_run_parameters(run_parameters)

SELECT = {
    "cluster_nmf": nmf,
    "cc_cluster_nmf": cc_nmf,
    "net_cluster_nmf": net_nmf,
    "cc_net_cluster_nmf": cc_net_nmf}


def main():
    import sys
    from knpackage.toolbox import get_run_directory_and_file
    from knpackage.toolbox import get_run_parameters

    run_directory, run_file = get_run_directory_and_file(sys.argv)
    run_parameters = get_run_parameters(run_directory, run_file)
    SELECT[run_parameters["method"]](run_parameters)


if __name__ == "__main__":
    main()
