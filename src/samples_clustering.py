"""
Created on Wed Jul 20 14:47:45 2016
@author: The KnowEnG dev team
"""

def nmf(run_parameters):
    '''nmf clustering'''
    from samples_clustering_toolbox import run_nmf
    run_nmf(run_parameters) 

def cc_nmf(run_parameters):
    '''kmeans consensus clustering of the nmf-based clusters'''
    from samples_clustering_toolbox import run_cc_nmf
    run_cc_nmf(run_parameters)

def net_nmf(run_parameters):
    '''net-nmf clustering "'''
    from samples_clustering_toolbox import run_net_nmf
    run_net_nmf(run_parameters)

def cc_net_nmf(run_parameters):
    '''kmeans consensus clustering of the net-nmf-based clusters'''
    from samples_clustering_toolbox import run_cc_net_nmf
    run_cc_net_nmf(run_parameters)

SELECT = { "nmf"       :nmf
         , "cc_nmf"    :cc_nmf
         , "net_nmf"   :net_nmf
         , "cc_net_nmf":cc_net_nmf }

def main():
    """
    This is the main function to perform sample clustering
    """
    import sys
    from knpackage.toolbox import get_run_directory_and_file
    from knpackage.toolbox import get_run_parameters
    
    run_directory, run_file = get_run_directory_and_file(sys.argv)
    run_parameters          = get_run_parameters(run_directory, run_file)

    SELECT[run_parameters["method"]](run_parameters)

if __name__ == "__main__":
    main()
