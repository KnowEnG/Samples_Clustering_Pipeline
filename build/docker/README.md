# Building The Samples Clustering Pipeline Docker Image

The Dockefile in this directory contains all the commands, in order, needed to build the **Samples Clustering Pipeline** docker image.

* Run the "make" command to build the **Samples Clustering Pipeline** docker image (output: docker image called "samples_clustering_pipeline" and a tag with today's date and time):
```
    make build_docker_image
```

* Login to docker hub. When prompted, enter your password and press enter:
```
    make login_to_dockerhub username=*enter your docker login here* email=*enter your email here*
```

* Upload your image to docker hub:
```
    make push_to_dockerhub
```

* * * 
## How to run this docker image
* * * 

### * Change directory to the directory  where you want to run

### * Run the following docker command to run the specific docker image:
```
docker run -v \`pwd\`:/home/test/run_dir/ -it knowengdev/samples_clustering_pipeline:10_18_2016 
```

### * Change directory to the test
```
cd test
```
### * Create a local directory "run_dir" and place all the run files in it
```
make env_setup
```

### * Run the Sample Clustering Pipeline
```
make run_cc_net_nmf
```
