# Building The Samples Clustering Pipeline Docker Image
The Dockefile in this directory contains all the commands, in order, needed to build samples_clustering_pipeline docker image.

Simply run the "make" command to build the samples_clustering_pipeline docker image. The makefile assumes there is a Dockerfile in current directory. 
```
    make build_docker_image
```
The results of the "make" command are a docker image called "samples_clustering_pipeline" and a tag with today's date and time.


Then login to docker hub before you push to it. When prompted, enter your password and press enter.
```
    make login_to_docker username=your docker login here email=your email as in your Docker profile here
```
Last, upload your image to docker hub!
```
    make push_to_dockerhub
```

