MKDIR_P = mkdir -p
RMR = rm -rf
GZIP = gzip
GZIP_D = gzip -d
INPUT_DATA_DIR = ./input_data
RUN_DIR = ./run_dir
RESULTS_DIR = $(RUN_DIR)/results
BENCHMARKS = ./test/benchmarks/
SCRIPT = ./src/samples_clustering.py


all: decompress_input_data create_run_dir copy_run_files

run_nmf:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file cluster_nmf_run_file.yml 

run_net_nmf:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file net_cluster_nmf_run_file.yml

run_cc_nmf:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file cc_cluster_nmf_run_file.yml 

run_cc_net_nmf:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file cc_net_cluster_nmf_run_file.yml


decompress_input_data:
	$(GZIP_D) $(INPUT_DATA_DIR)/*
 
compress_input_data:
	$(GZIP) $(INPUT_DATA_DIR)/*

create_input_data:
	$(MKDIR_P) $(INPUT_DATA_DIR)	

create_run_dir:
	$(MKDIR_P) $(RESULTS_DIR) 

copy_run_files:
	cp $(BENCHMARKS)/*.yml $(RUN_DIR) 

clean_dir_recursively:
	$(RMR) $(RUN_DIR)