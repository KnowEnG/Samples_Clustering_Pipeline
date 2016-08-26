# Samples_Clustering_Pipeline
This pipeline selects one of three methods to **rank** a user supplied gene set **against** a KnowEnG gene sets collection

## Steps to run pipelines
###1. Setup github access:
__Access__ KnowEnG-Research github repo

###2. Get a copy of the Samples_Clustering_Pipeline code, data
__Run__ the following command to get Samples_Clustering_Pipeline repo
```
 git clone https://github.com/KnowEnG-Research/Samples_Clustering_Pipeline.git
```
 
###3. Configure your environment to have the following packages
  ```
 System: ubuntu:14.04
 apt-get install -y python3-pip
 apt-get install -y libblas-dev liblapack-dev libatlas-base-dev gfortran
 pip3 install -I numpy==1.11.1
 pip3 install -I pandas==0.18.1 
 pip3 install -I scipy==0.18.0
 pip3 install -I scikit-learn==0.17.1
 pip3 install -I vcversioner==2.16.0.0
 pip3 install -I knpackage==0.1.2
 apt-get install -y libfreetype6-dev libxft-dev 
 pip3 install -I matplotlib==1.4.2
 pip3 install pyyaml
```

###4. start in the pipeline repo home directory

```
cd Samples_Clustering_Pipeline
```

 
###5. Run makefile targets
  * Prepare input data and running directories. 
 ```
  make preparation
 ```
 
  * Run the pipeline you desire
 ```
make run_nmf
make run_net_nmf
make run_cc_nmf
make run_cc_net_nmf
 ```
 
  * Clean the running environment and revert the local repository to original state after you finish your tests and analysis
 ```
  make final_clean 
 ```
 

###6. Run methods seperately

* Create your own run directory outside Samples_Clustering_Pipeline repo
 ```
 mkdir run_dir
 ```

* Create results directory to save output files under run directory
 ```
 cd run_dir
 mkdir results
 ```
 
####Make sure you are in the run_dir directory.

### non-negative matrix factorization (nmf)
1. Copy `cluster_nmf_run_file.yml` into run_dir
  ```
  cp ../Samples_Clustering_Pipeline/test/benchmarks/cluster_nmf_run_file.yml cluster_nmf_run_file.yml
  ```
  
2. Make sure the directories of the input data in `cluster_nmf_run_file.yml` are correct
  
  pg_network_file_name:
  ```
  /../Samples_Clustering_Pipeline/input_data/final_clean_4col.edge
  ```
  samples_file_name:
  ```
  /../Samples_Clustering_Pipeline/input_data/final_clean_full_matrix.df
  ```
  
3. Run nmf
  ```
  export PYTHONPATH='../Samples_Clustering_Pipeline/src':$PYTHONPATH    
  
  python3 ../Samples_Clustering_Pipeline/src/samples_clustering.py -run_directory ./ -run_file cluster_nmf_run_file.yml
  ```
  
4. Output files are saved in results directory
  Generate `labels_data.tsv` file with timestamp.
