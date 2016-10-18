# Running the docker image to set up an environment

## Set up and run in a terminal (if you have docker installed):
1 Change directory to the directory  where you want to run.

2 docker run -v \`pwd\`:\`pwd\` -it knowengdev/samples_clustering_pipeline:09_01_2016

3 make all

4 make run_cc_net_nmf

* The make options in Samples_Clustering_Pipeline/README.md apply.

* Check on docker.hub to get the latest image. 

* If you don't "cp" your data into the volume you mounted it will disappear when you exit docker.

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
