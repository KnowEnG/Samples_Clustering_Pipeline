# ----------------------------------------------------------------
# Variables setup Section                                        -
# ----------------------------------------------------------------
SCRIPT = 		../src/samples_clustering.py
RUN_DIR = 		./run_dir
RUN_FILES = 		../data/run_files
RESULTS_DIR = 		$(RUN_DIR)/results
INPUT_DATA_DIR = 	../data
UNIT_TEST_TEMP_DIR =	 ./unit/tmp

RMR = rm -rf
GZIP = gzip
GZIP_D = gzip -d
MKDIR_P = mkdir -p

# ----------------------------------------------------------------
# Benchmarks setup Section                                       -
# ----------------------------------------------------------------
env_setup: create_run_dir copy_run_files

decompress_input_data:
	$(GZIP_D) $(INPUT_DATA_DIR)/networks/* $(INPUT_DATA_DIR)/spreadsheets/*
 
compress_input_data:
	$(GZIP) $(INPUT_DATA_DIR)/networks/* $(INPUT_DATA_DIR)/spreadsheets/*

create_run_dir:
	$(MKDIR_P) $(RESULTS_DIR) 

copy_run_files:
	cp $(RUN_FILES)/*.yml $(RUN_DIR) 

clean_dir_recursively:
	$(RMR) $(RUN_DIR)

# ----------------------------------------------------------------
# Large Benchmarks RUN Section                                   -
# ----------------------------------------------------------------
run_nmf:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file nmf.yml 

run_net_nmf:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file net_nmf.yml

run_cc_nmf:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file cc_nmf.yml 

run_cc_net_nmf:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file cc_net_nmf.yml

run_all_methods: run_nmf run_net_nmf run_cc_nmf run_cc_net_nmf

# ----------------------------------------------------------------
# - UNIT TESTS RUN SECTION                                       -
# ----------------------------------------------------------------
unit_test_setup:
	$(MKDIR_P) $(UNIT_TEST_TEMP_DIR)

unit_tests:
	cd unit; make all_unit_tests
	
run_benchmark_test:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file cc_net_nmf_paralell.yml	

# ----------------------------------------------------------------
# - CLEAN UP SECTION                                             -
# ----------------------------------------------------------------	
unit_test_cleanup:
	$(RMR) $(UNIT_TEST_TEMP_DIR)


