### Answers data UCSD_precalc_UCEC_ST90.tsv
* The precalculated data is in: data/verification/UCSD_precalc_UCEC_ST90.tsv
* The run_benchmark_test target in the ./test/Makefile runs the benchmark_cc_net_nmf_run_file.yml
* The members of each cluster may have different numbering assignments (0,2,1) vs (1,2,3).
* The parameter "number_of_bootstraps" should be at least 200, 1000 would be better
* The parameter "processing_method" must be set for your environment: dist_comp is fastest, parl_loc is fast, and serial works anywhere
