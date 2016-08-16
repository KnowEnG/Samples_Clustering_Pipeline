# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 14:47:45 2016

@author: The Gene Sets Characterization dev team

"""
import time
import os
def disp_run_parameters(run_parameters):
    if int(run_parameters['verbose']) == 0:
        return
    else:
        nstr = time.strftime("%a_%d_%b_%Y_%H_%M_%S", time.localtime())
        print('\n\nBegin {}\t{}'.format(run_parameters["method"], nstr))
        for p in run_parameters:
            print('{} : {}'.format(p, run_parameters[p]))
        print('\n')

def log_run_parameters(run_parameters):
    nstr = time.strftime("%a_%d_%b_%Y_%H_%M_%S", time.localtime())
    full_file = os.path.join(run_parameters["run_directory"], run_parameters["method"] + nstr + '.log')
    run_parameters[run_parameters["method"]] = nstr
    with open(full_file, 'a') as log:
        for p in run_parameters:
            log.write('{} : {}\n'.format(p, run_parameters[p]))

def nmf(run_parameters):
    '''nmf clustering'''
    from sample_clustering_toolbox import run_nmf
    run_parameters['display_clusters'] = 0
    disp_run_parameters(run_parameters)
    t0 = time.time()
    run_nmf(run_parameters)
    run_parameters['running_time'] = time.time() - t0
    log_run_parameters(run_parameters)
    print('\nrun_nmf total time: {}\n'.format(time.time() - t0))

def cc_nmf(run_parameters):
    '''kmeans consensus clustering of the nmf-based clusters'''
    from sample_clustering_toolbox import run_cc_nmf
    run_parameters['display_clusters'] = 0
    disp_run_parameters(run_parameters)
    t0 = time.time()
    run_cc_nmf(run_parameters)
    run_parameters['running_time'] = time.time() - t0
    log_run_parameters(run_parameters)
    print('\nrun_cc_nmf total time: {}\n'.format(time.time() - t0))

def net_nmf(run_parameters):
    '''net-nmf clustering "'''
    from sample_clustering_toolbox import run_net_nmf
    run_parameters['display_clusters'] = 0
    disp_run_parameters(run_parameters)
    t0 = time.time()
    run_net_nmf(run_parameters)
    run_parameters['running_time'] = time.time() - t0
    log_run_parameters(run_parameters)
    print('\nrun_net_nmf total time: {}\n'.format(time.time() - t0))

def cc_net_nmf(run_parameters):
    '''kmeans consensus clustering of the net-nmf-based clusters'''
    from sample_clustering_toolbox import run_cc_net_nmf
    run_parameters['display_clusters'] = 0
    disp_run_parameters(run_parameters)
    t0 = time.time()
    run_cc_net_nmf(run_parameters)
    run_parameters['running_time'] = time.time() - t0
    log_run_parameters(run_parameters)
    print('\nrun_nmf total time: {}\n'.format(time.time() - t0))

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
