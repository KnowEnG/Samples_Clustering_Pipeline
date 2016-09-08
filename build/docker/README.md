# Samples Clustering Pipeline
This Dockefile contains all the commands, in order, needed to build Samples Clustering Pipeline docker image.

## Getting Started
Simply run the "make" command to build the samples_clustering_pipeline image. The makefile assumes there is a Dockerfile in current directory. The results of the "make" command are a docker image called "sample_clustering_pipeline" and a tag with today's date and time.
```
    make build_docker_image
```
Then login to docker hub before you push to it. When prompted, enter your password and press enter.
```
    make login_to_docker username=your docker login here email=your email as in your Docker profile here
```
Last, upload your image to docker hub!
```
    make push_to_dockerhub
```

