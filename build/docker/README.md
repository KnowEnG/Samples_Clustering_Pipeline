# Samples Clustering Pipeline
This Dockefile contains all the commands, in order, needed to build Samples Clustering Pipeline docker image.

## Getting Started
Simply run the following command to build an image with the latest code changes. It uses the Dockerfile in current directory 
and generates a docker image named sample_clustering_pipelin with the tag indicating the date what is was created.
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

