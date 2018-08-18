"""
sobh@illinois.edu

"""

import filecmp
import os
import time

verification_dir = '../data/verification/'
results_dir = '../test/run_dir/results'


def verify_benchmark(option, BENCHMARK_name_list, BENCHMARK_YML):
    run_command = 'python3 ../src/samples_clustering.py -run_directory ./run_dir -run_file ' + BENCHMARK_YML
    os.system(run_command)

    All_files_in_results_dir = os.listdir(results_dir)

    num_failed_tests = 0
    num_succeed_tests = 0
    for f in All_files_in_results_dir:
        for BENCHMARK_name in BENCHMARK_name_list:
            if BENCHMARK_name in f:
                RESULT    = os.path.join(results_dir, f)
                BENCHMARK = os.path.join(verification_dir, option, BENCHMARK_name + '.tsv')
                if filecmp.cmp(RESULT, BENCHMARK) == True:
                    num_succeed_tests += 1
                    print(BENCHMARK, '______ PASS ______')
                else:
                    num_failed_tests += 1
                    print(BENCHMARK, '****** FAIL ******')
    return num_succeed_tests, num_failed_tests


def main():
    BENCHMARK = {
        'nmf': ['BENCHMARK_1_SC_nmf.yml',
                'clustering_evaluation_result_nmf',
                'consensus_matrix_nmf',
                'genes_averages_by_cluster_nmf',
                'genes_by_samples_heatmap_nmf',
                'genes_variance_nmf',
                'samples_label_by_cluster_nmf',
                'silhouette_overall_score_nmf',
                'silhouette_per_cluster_score_nmf',
                'silhouette_per_sample_score_nmf',
                'top_genes_by_cluster_nmf'
        ],
        'net_nmf': [
                'BENCHMARK_2_SC_net_nmf.yml',
                'clustering_evaluation_result_net_nmf',
                'consensus_matrix_net_nmf',
                'genes_averages_by_cluster_net_nmf',
                'genes_by_samples_heatmap_net_nmf',
                'genes_variance_net_nmf',
                'samples_label_by_cluster_net_nmf',
                'silhouette_overall_score_net_nmf',
                'silhouette_per_cluster_score_net_nmf',
                'silhouette_per_sample_score_net_nmf',
                'top_genes_by_cluster_net_nmf'
        ],
        'cc_nmf': [
                'BENCHMARK_4_SC_cc_nmf_parallel_shared.yml',
                'clustering_evaluation_result_cc_nmf',
                'consensus_matrix_cc_nmf',
                'genes_averages_by_cluster_cc_nmf',
                'genes_by_samples_heatmap_cc_nmf',
                'genes_variance_cc_nmf',
                'samples_label_by_cluster_cc_nmf',
                'silhouette_overall_score_cc_nmf',
                'silhouette_per_cluster_score_cc_nmf',
                'silhouette_per_sample_score_cc_nmf',
                'top_genes_by_cluster_cc_nmf'
        ],
        'cc_net_nmf': [
                'BENCHMARK_7_SC_cc_net_nmf_parallel_shared.yml',
                'clustering_evaluation_result_cc_net_nmf',
                'consensus_matrix_cc_net_nmf',
                'genes_averages_by_cluster_cc_net_nmf',
                'genes_by_samples_heatmap_cc_net_nmf',
                'genes_variance_cc_net_nmf',
                'samples_label_by_cluster_cc_net_nmf',
                'silhouette_overall_score_cc_net_nmf',
                'silhouette_per_cluster_score_cc_net_nmf',
                'silhouette_per_sample_score_cc_net_nmf',
                'top_genes_by_cluster_cc_net_nmf']
    }
    os.system('make env_setup')
    start_time = time.time()

    total_success = 0
    total_failure = 0
    for option in BENCHMARK.keys():

        BENCHMARK_list = BENCHMARK[option]
        BENCHMARK_YML  = BENCHMARK_list[0]
        print()
        print("INFO: Running test ", "./run_dir/results/" + BENCHMARK_YML)

        # for BENCHMARK_name in BENCHMARK_list[1:]:

        num_succeed_tests, num_failed_tests = verify_benchmark(option, BENCHMARK_list[1:], BENCHMARK_YML)

        total_success += num_succeed_tests
        total_failure += num_failed_tests

        os.system('rm ./run_dir/results/*')

    end_time = time.time()

    print()
    print("Ran {} tests in {}s".format(total_success + total_failure, end_time - start_time))
    if (total_failure == 0):
        print("OK")
        print()
    else:
        print("FAILED(errors={})".format(total_failure))


if __name__ == "__main__":
    main()
