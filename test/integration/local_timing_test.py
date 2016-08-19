
"""
    script for time testing the clustering mehtods in sample_clustering_toolbox.py
    lanier4@illinois.edu
    or the KnowEnG clustering team
"""
import os
import time
from knpackage.toolbox import get_run_parameters
from sample_clustering_toolbox import run_cc_net_nmf
from sample_clustering_toolbox import run_cc_nmf
from sample_clustering_toolbox import run_net_nmf
from sample_clustering_toolbox import run_nmf

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

local_run_dir_path = '/Users/lanier4/BigDataTank/nbs_run_silent'
cc_net_nmf_file = 'cc_net_cluster_nmf_run_file.yml'
cc_nmf_file = 'cc_cluster_nmf_run_file.yml'
net_nmf_file = 'net_cluster_nmf_run_file.yml'
nmf_file = 'cluster_nmf_run_file.yml'

t_overall_start = time.time()

# run consensus clustering networked nmf
run_filename = cc_net_nmf_file

run_parameters = get_run_parameters(local_run_dir_path, run_filename)
run_parameters['verbose'] = 1
run_parameters['display_clusters'] = 0
disp_run_parameters(run_parameters)
t0 = time.time()
run_cc_net_nmf(run_parameters)
run_parameters['running_time'] = time.time() - t0
log_run_parameters(run_parameters)
print('\n\t{} complete in {} seconds'.format(run_filename, run_parameters['running_time']))


# run consensus clustering nmf
run_filename = cc_nmf_file

run_parameters = get_run_parameters(local_run_dir_path, run_filename)
run_parameters['verbose'] = 1
run_parameters['display_clusters'] = 0
disp_run_parameters(run_parameters)
t0 = time.time()
run_cc_nmf(run_parameters)
run_parameters['running_time'] = time.time() - t0
log_run_parameters(run_parameters)
print('\n\t{} complete in {} seconds'.format(run_filename, run_parameters['running_time']))

# run consensus clustering networked nmf
run_filename = net_nmf_file

run_parameters = get_run_parameters(local_run_dir_path, run_filename)
run_parameters['verbose'] = 1
run_parameters['display_clusters'] = 0
disp_run_parameters(run_parameters)
t0 = time.time()
run_net_nmf(run_parameters)
run_parameters['running_time'] = time.time() - t0
log_run_parameters(run_parameters)
print('\n\t{} complete in {} seconds'.format(run_filename, run_parameters['running_time']))

# run consensus clustering networked nmf
run_filename = nmf_file

run_parameters = get_run_parameters(local_run_dir_path, run_filename)
run_parameters['verbose'] = 1
run_parameters['display_clusters'] = 0
disp_run_parameters(run_parameters)
t0 = time.time()
run_nmf(run_parameters)
run_parameters['running_time'] = time.time() - t0
log_run_parameters(run_parameters)
print('\n\t{} complete in {} seconds'.format(run_filename, run_parameters['running_time']))

print('\n\tAll Tests complete in {} seconds'.format(time.time() - t_overall_start))
