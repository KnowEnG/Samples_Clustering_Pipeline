# KnowEnG's Samples Clustering Pipeline 
This is the Knowledge Engine for Genomics (KnowEnG), an NIH BD2K Center of Excellence, Samples Clustering Pipeline.

This pipeline **clusters** the columns of a given spreadsheet, where spreadsheet's columns correspond to sample-labels and rows correspond to gene-labels.

There are four clustering methods that one can choose from:


| **Options**                                      | **Method**                           | **Parameters** |
| ------------------------------------------------ | -------------------------------------| -------------- |
| Clustering                                       | nmf                                  | nmf            |
| Consensus Clustering                             | bootstrapping with nmf               | cc_nmf         |
| Clustering with network regularization           | network-based nmf                    | net_nmf        |
| Consensus Clustering with network regularization | bootstrapping with network-based nmf | cc_net_nmf     |


Note: all of the clustering methods mentioned above use the non-negative matrix factorization (nmf) as the main clustering algorithm.


If a pheotype data file is include pipeline **evaluates** the clustering result of KnowEnG's Samples Clustering Pipeline.

There are two evaluation methods:

| **Method**                                      | **Trait Type**                          |
| ------------------------------------------------ | ------------------------------------- |
| one-way ANOVA(f_oneway)                               | Continuous                                | 
| one-way chi square test(chisquare)                                     | Categorical          |

* * * 
## How to run this pipeline with Our data
* * * 
###1. Clone the Samples_Clustering_Pipeline Repo
```
 git clone https://github.com/KnowEnG-Research/Samples_Clustering_Pipeline.git
```
 
###2. Install the following (Ubuntu or Linux)
  ```
 pip3 install pyyaml
 pip3 install knpackage
 pip3 install scipy==0.18.0
 pip3 install numpy==1.11.1
 pip3 install pandas==0.18.1
 pip3 install matplotlib==1.4.2
 pip3 install scikit-learn==0.17.1
 
 apt-get install -y python3-pip
 apt-get install -y libfreetype6-dev libxft-dev
 apt-get install -y libblas-dev liblapack-dev libatlas-base-dev gfortran
```

###3. Change directory to Samples_Clustering_Pipeline

```
cd Samples_Clustering_Pipeline
```

###4. Change directory to test

```
cd test
```
 
###5. Create a local directory "run_dir" and place all the run files in it
```
make env_setup
```

###6. Use one of the following "make" commands to select and run a clustering option:


| **Command**         | **Option**                                       | 
|:------------------- |:------------------------------------------------ | 
| make run_nmf        | Clustering                                       |
| make run_net_nmf     | Clustering with network regularization           |
| make run_cc_nmf_serial     | Consensus Clustering                             |
| make run_cc_nmf_parallel_shared     | Consensus Clustering                             |
| make run_cc_net_nmf_serial | Consensus Clustering with network regularization |
| make run_cc_net_nmf_parallel_shared | Consensus Clustering with network regularization |

 
* * * 
## How to run this pipeline with Your data
* * * 

__***Follow steps 1-3 above then do the following:***__

### * Create your run directory

 ```
 mkdir run_directory
 ```

### * Change directory to the run_directory

 ```
 cd run_directory
 ```

### * Create your results directory

 ```
 mkdir results_directory
 ```
 
### * Create run_paramters file  (YAML Format)
 ``` 
 Look for examples of run_parameters in the Sample_Clustering_Pipeline/data/run_files TEMPLATE_cc_net_cluster_nmf.yml
 ```
### * Modify run_paramters file  (YAML Format)
Change processing_method to one of: serial, parallel depending on your machine.
```
processing_method: serial
```

set the data file targets to the files you want to run, and the parameters as appropriate for your data.


### * Run the Samples Clustering Pipeline:

  * Update PYTHONPATH enviroment variable
   ``` 
   export PYTHONPATH='../src':$PYTHONPATH    
   ```
   
  * Run
   ```
  python3 ../src/samples_clustering.py -run_directory ./ -run_file TEMPLATE_cc_net_cluster_nmf.yml
   ```

* * * 
## Description of "run_parameters" file
* * * 

| **Key**                   | **Value** | **Comments** |
| ------------------------- | --------- | ------------ |
| method                    | **nmf**, **cc_nmf**, **net_nmf** or **cc_net_nmf**  | Choose clustering method |
| gg_network_name_full_path | directory+gg_network_name |Path and file name of the 4 col network file |
| spreadsheet_name_full_path | directory+spreadsheet_name|  Path and file name of user supplied gene sets |
| phenotype_data_full_path | directory+phenotype_data_name| Path and file name of user supplied phenotype data |
| threshold | 10 | cluster eval - catagorical vs continuous cut off level |
| results_directory | directory | Directory to save the output files |
| tmp_directory | directory | Directory to save the intermediate files |
| rwr_max_iterations | 100| Maximum number of iterations without convergence in random walk with restart |
| rwr_convergence_tolerence | 1.0e-8 | Frobenius norm tolerence of spreadsheet vector in random walk|
| rwr_restart_probability | 0.7 | alpha in `V_(n+1) = alpha * N * Vn + (1-alpha) * Vo` |
| rows_sampling_fraction| 0.8| Select 80% of spreadsheet rows|
| cols_sampling_fraction| 0.8| Select 80% of spreadsheet columns|
| number_of_bootstraps| 4 | Number of random samplings |
| number_of_clusters| 3 | Estimated number of clusters |
| nmf_conv_check_freq| 50 | Check convergence at given frequency |
| nmf_max_invariance| 200 | Maximum number of invariance |
| nmf_max_iterations| 10000 | Maximum number of iterations |
| nmf_penalty_parameter| 1400 | Penalty parameter |
| top_number_of_genes| 100 | Number of top genes selected |
| processing_method| serial or parallel or distribute | Choose processing method |
| parallelism| number of cores to use in parallel processing | Set number of cores for speed or memory |

gg_network_name = STRING_experimental_gene_gene.edge</br>
spreadsheet_name = ProGENI_rwr20_STExp_GDSC_500.rname.gxc.tsv</br>
phenotype_data_name = UCEC_phenotype.txt

* * * 
## Description of Output files saved in results directory
* * * 

* Output files of all four methods save genes by sample heatmap variances per row with name **genes_variance_{method}_{timestamp}_viz.tsv**.</br>

 |  |**variance**|
 | :--------------------: |:--------------------:|
 | **gene 1**|float|
 |...|...|
 | **gene m**| float|

* Output files of all four methods save genes by samples heatmap with name **genes_by_samples_heatmp_{method}_{timestamp}_viz.tsv**.</br>

 |  |**sample 1**|...|**sample n**|
 | :--------------------: |:--------------------:|:--------------------:|:--------------------:|
 | **gene 1**|float|...|float|
 |...|...|...|...|
 | **gene m**|float|...|float|

* Output files of all four methods save samples by samples heatmap with name **consensus_matrix_{method}_{timestamp}_viz.tsv**.</br>

 |  |**sample 1**|...|**sample n**|
 | :--------------------: |:--------------------:|:--------------------:|:--------------------:|
 | **sample 1**|float|...|float|
 |...|...|...|...|
 | **sample n**|float|...|float|
 
* Output files of all four methods save patients to cluster map with name **samples_labeled_by_cluster_{method}_{timestamp}_viz.tsv**.</br>

 |    |**cluster**|
 | :--------------------: |:--------------------:|
 | **sample 1** | int|
 |...|...|
 | **sample n** |int|
 
* Output files of all four methods save gene scores by cluster with name **genes_averages_by_cluster_{method}_{timestamp}_viz.tsv**.</br>

 |  |**cluster 1**|...|**cluster k**|
 | :--------------------: |:--------------------:|:--------------------:|:--------------------:|
 | **gene 1**|float|...|float|
 |...|...|...|...|
 | **gene m**|float|...|float|
 
* Output files of all four methods save spreadsheet with top ranked genes per sample with name **top_genes_by_cluster_{method}_{timestamp}_download.tsv**.</br>

 |  |**cluster 1**|...|**cluster k**|
 | :--------------------: |:--------------------:|:--------------------:|:--------------------:|
 | **gene 1**|1/0|...|1/0|
 |...|...|...|...|
 | **gene m**|1/0|...|1/0|
  
* All four methods save **silhouette number of clusters** and **corresponding silhouette score** with name silhouette_average\_{method}\_{timestamp}\_viz.tsv.</br>
 ```
 File Example: 
 silhouette number of clusters = 3, corresponding silhouette score = 1
 ```

* Output files of all four methods save patients to cluster map with name **phenotypes_labeled_by_cluster_{method}_{timestamp}_viz.tsv**.</br>

 | **sample id** |**cluster**|**phenotype 1**|...|**phenotype k**|
 | :--------------------: |:--------------------:|:--------------------:|:--------------------:|:--------------------:|
 | **sample 1**|int|mixed type|...|mixed type|
 |...|...|...|...|...|
 | **sample n**|int|mixed type|...|mixed type|
 
 
 * The clustering evaluation output file has the name 
 **clustering_evaluation_result_{timestamp}.tsv**.</br>

 |  |**Measure**|**Trait_length_after_dropna**| **Sample_number_after_dropna**|**chi/fval**|**pval**|
 | :--------------------: |:--------------------:|:--------------------:|:--------:|:-------:|:--------------------:|
 | **sample 1**|f_oneway|int(more than threshold)|int|float|float|
 |...|...|...|...|...|...|
 | **sample m**|chisquare|int(less than threshold)|int|float|float|
