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

# Running the docker image to set up an environment

## Set up and run in a terminal (if you have docker installed):
1 Change directory to the directory  where you want to run.

2 docker run -v \`pwd\`:/home/test/run_dir/ -it knowengdev/samples_clustering_pipeline:10_18_2016 

3 cd test

4 make env_setup

5 edit the .yml file (use the comments to see options)

6 make run_cc_net_nmf

* The make options in Samples_Clustering_Pipeline/README.md apply.

* Check on docker.hub to get the latest image. 

* If you don't "cp" your data into the volume you mounted it will disappear when you exit docker.
