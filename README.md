# KnowEnG's Samples Clustering Pipeline 

*KnowEnG stands for Knowledge Engine for Genomics. KnowEnG is one of the NIH Big Data to Knowledge Centers of Excellence*

This pipeline clusters the columns of a given spreadsheet, where spreadsheet's columns correspond to sample-labels and rows correspond to gene-labels).

There are four possible clustering methods that one can chose from:


| **Options**                                      | **Method**                           | **Parameters** |
| ------------------------------------------------ | -------------------------------------| -------------- |
| Clustering                                       | nmf                                  | nmf            |
| Consensus Clustering                             | bootstrapping with nmf               | cc_nmf         |
| Clustering with network regularization           | network-based nmf                    | net_nmf        |
| Consensus Clustering with network regularization | bootstrapping with network-based nmf | cc_net_nmf     |


Note: all of the clustering methods mentioned above use the non-negative matrix factorization (nmf) as the main clustering algorithm.

## How to run this pipeline with provided data
###1. Get Access to KnowEnG-Research Repo:
Email omarsobh@illinois.edu infrastructure team (IST) lead to:

1. __Access__ KnowEnG-Research github repo

###2. Clone the Samples_Clustering_Pipeline Repo
```
 git clone https://github.com/KnowEnG-Research/Samples_Clustering_Pipeline.git
```
 
###3. Install the following (Mac OS or Linux)
  ```
 apt-get install -y python3-pip
 apt-get install -y libblas-dev liblapack-dev libatlas-base-dev gfortran
 pip3 install numpy==1.11.1
 pip3 install pandas==0.18.1
 pip3 install scipy==0.18.0
 pip3 install scikit-learn==0.17.1
 apt-get install -y libfreetype6-dev libxft-dev
 pip3 install matplotlib==1.4.2
 pip3 install pyyaml
 pip3 install knpackage
```

###4. Change directory to  the Samples_Clustering_Pipeline/test

```
cd Samples_Clustering_Pipeline/test
```

 
###5. Use the following "make" command to create a local directory "run_dir" and place all the run files in it.
 ```
  make env_setup
 ```

###6. Use one of the following "make" commands to select and run a clustering option:


| **Command**         | **Option**                                       | 
|:------------------- |:------------------------------------------------ | 
| make run_nmf        | Clustering                                       |
| make run_cc_nmf     | Consensus Clustering                             |
| make run_cc_nmf     | Clustering with network regularization           |
| make run_cc_net_nmf | Consensus Clustering with network regularization |

 
## How to run it with your data 
### Setup your run environment

* Create a  run directory

 ```
 mkdir run_directory_name
 ```

* Create results directory to save output files under run directory

 ```
 cd run_directory_name
 mkdir results_directory_name
 ```
 
* Create run_paramerters file (commented template file: data/run_files/zzz_clustering_run_file_template.yml) 

  custom_run_file.yml

* Make sure the directories of the input data in `custom_run_file.yml` are correct
 
* Run Samples Clustering Pipeline

```
  export PYTHONPATH='../Samples_Clustering_Pipeline/src':$PYTHONPATH    
  python3 ../Samples_Clustering_Pipeline/src/samples_clustering.py -run_directory ./ -run_file custom_run_file.yml
  ```
  
### Output files are saved in results_directory_name/ (CAUTION: make clean_dir_recursively will erase the ouput files)
 
## Note: docker image file instructions are in the build directory
