<<<<<<< HEAD
# Samples Clustering Pipeline 
This is the Knowledge Engine for Genomics (KnowEnG), an NIH BD2K Center of Excellence, samples clustering pipeline. 
This pipeline clusters samples in a user submitted spreadsheet (with samples as columns and genes as rows). 

One can select one of four clustering options that are based non-negative matrix factorization (nmf):


| **Options**                                      | **Method**                           | **Parameters** |
| ------------------------------------------------ | -------------------------------------| -------------- |
| Clustering                                       | nmf                                  | nmf            |
| Consensus Clustering                             | bootstrapping with nmf               | cc_nmf         |
| Clustering with network regularization           | network-based nmf                    | net_nmf        |
| Consensus Clustering with network regularization | bootstrapping with network-based nmf | cc_net_nmf     |

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

###4. Change directory to  the Samples_Clustering_Pipeline

```
cd Samples_Clustering_Pipeline
```

 
###5. Use the following "make" command to create a local directory "run_dir" and place all the parameters files in it.
 ```
  make run_env
 ```

###6. Use one of the following "make" commands to select and run a clustering option:


| **Command**         | **Option**                                       | 
| ------------------- | -------------------------------------------------| 
| make run_nmf        | Clustering                                       |
| make run_cc_nmf     |  Consensus Clustering                            |
| make run_cc_nmf     | Clustering with network regularization           |
| make run_cc_net_nmf | Consensus Clustering with network regularization |

 
## How to run it with your data 
###6. Setup your run environment

* Create a  run directory

 ```
 mkdir run_directory_name
 ```

* Create results directory to save output files under run directory

 ```
 cd run_directory_name
 mkdir results_directory_name
 ```
 
* Create run_paramerters file (yml format) 

  filename.yml

* Make sure the directories of the input data in `cluster_nmf_run_file.yml` are correct
 
* Run Samples Clustering Pipeline

```
  export PYTHONPATH='../Samples_Clustering_Pipeline/src':$PYTHONPATH    
  python3 ../Samples_Clustering_Pipeline/src/samples_clustering.py -run_directory ./ -run_file file_name.yml
  ```
  
* Output files are saved in results directory
 
## How to set up and run in a terminal if you have docker installed:
1 Change directory to the directory where you want to run.

2 docker run -v `pwd`:`pwd` -it knowengdev/samples_clustering_pipeline:09_01_2016

3 make all

4 make run_cc_net_nmf

* The make instructions above apply.

* Check on docker.hub to get the latest image. 

* If you don't "cp" your data into the volume you mounted it will disappear when you exit docker.
=======
# Samples_Clustering_Pipeline
Clusters a given sample using one of four clustering options

* User submits a spreadsheet with samples as columns and genes as rows. 
* System clusters samples.
* Specific example: network-based stratification of patients 
>>>>>>> upstream/master
