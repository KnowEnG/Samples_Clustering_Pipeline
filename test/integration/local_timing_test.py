
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

WRITE_LOGS = True
LOGFILE_NAME = 'local_timing'
DISPLAY_PARS = True

def disp_closing_time(run_file, run_time):
    print('\n\t\t^ Parameters Dictionary: {} run in {} sec.\n'.format(run_file, run_time))

def disp_run_parameters(run_parameters):
    if DISPLAY_PARS == 0:
        return
    else:
        nstr = time.strftime("%a_%d_%b_%Y_%H_%M_%S", time.localtime())
        print('\n\tBegin {}\t{}. Using Prarmeters:\n'.format(run_parameters["method"], nstr))
        keys = list(run_parameters.keys())
        keys.sort()
        for k in keys:
            print('{} : {}'.format(k, run_parameters[k]))
        print('\n')

def log_run_parameters(run_parameters):
    if WRITE_LOGS:
        full_file = run_parameters['LogFile']
        with open(full_file, 'a') as log:
            log.write('\n')
            keys = list(run_parameters.keys())
            keys.sort()
            for k in keys:
                log.write('{} : {}\n'.format(k, run_parameters[k]))
            if len(run_parameters) > 2:
                log.write('\n\t^ Section Running time was {}\n\n'.format(run_parameters['running_time']))
    else:
        pass

cc_net_nmf_file = 'cc_net_cluster_nmf_run_file.yml'
cc_nmf_file = 'cc_cluster_nmf_run_file.yml'
net_nmf_file = 'net_cluster_nmf_run_file.yml'
nmf_file = 'cluster_nmf_run_file.yml'

local_run_dir_path = '/Users/lanier4/BigDataTank/nbs_run_silent'
if WRITE_LOGS:
    nstr = time.strftime("%a_%d_%b_%Y_%H_%M_%S", time.localtime())
    full_log_file = os.path.join(local_run_dir_path, LOGFILE_NAME + nstr + '.log')
    print('\n\t\tLog file: {}\n'.format(full_log_file))
else:
    full_log_file = 'No Log File'
    print('\n\t\t{}\n'.format(full_log_file))

t_overall_start = time.time()

# run consensus clustering networked nmf
run_filename = cc_net_nmf_file

run_parameters = get_run_parameters(local_run_dir_path, run_filename)
run_parameters['Actual_Parameters_File'] = os.path.join(run_parameters['run_directory'], run_filename)
run_parameters['verbose'] = 1
run_parameters['display_clusters'] = 0
disp_run_parameters(run_parameters)
t0 = time.time()
run_cc_net_nmf(run_parameters)

run_parameters['running_time'] = time.time() - t0
run_parameters['LogFile'] = full_log_file
log_run_parameters(run_parameters)
disp_closing_time(run_filename, run_parameters['running_time'])

# run consensus clustering nmf
run_filename = cc_nmf_file

run_parameters = get_run_parameters(local_run_dir_path, run_filename)
run_parameters['Actual_Parameters_File'] = os.path.join(run_parameters['run_directory'], run_filename)
run_parameters['verbose'] = 1
run_parameters['display_clusters'] = 0
disp_run_parameters(run_parameters)
t0 = time.time()
run_cc_nmf(run_parameters)

run_parameters['running_time'] = time.time() - t0
run_parameters['LogFile'] = full_log_file
log_run_parameters(run_parameters)
disp_closing_time(run_filename, run_parameters['running_time'])

# run consensus clustering networked nmf
run_filename = net_nmf_file

run_parameters = get_run_parameters(local_run_dir_path, run_filename)
run_parameters['Actual_Parameters_File'] = os.path.join(run_parameters['run_directory'], run_filename)
run_parameters['verbose'] = 1
run_parameters['display_clusters'] = 0
disp_run_parameters(run_parameters)
t0 = time.time()
run_net_nmf(run_parameters)

run_parameters['running_time'] = time.time() - t0
run_parameters['LogFile'] = full_log_file
log_run_parameters(run_parameters)
disp_closing_time(run_filename, run_parameters['running_time'])

# run consensus clustering networked nmf
run_filename = nmf_file

run_parameters = get_run_parameters(local_run_dir_path, run_filename)
run_parameters['Actual_Parameters_File'] = os.path.join(run_parameters['run_directory'], run_filename)
run_parameters['verbose'] = 1
run_parameters['display_clusters'] = 0
disp_run_parameters(run_parameters)
t0 = time.time()
run_nmf(run_parameters)

run_parameters['running_time'] = time.time() - t0
run_parameters['LogFile'] = full_log_file
log_run_parameters(run_parameters)
disp_closing_time(run_filename, run_parameters['running_time'])


cap_time = time.time() - t_overall_start
run_parameters = {'LogFile':full_log_file,'running_time':cap_time}
log_run_parameters(run_parameters)
disp_closing_time('All Tests complete', cap_time)
