FROM knowengdev/base_image:09_01_2016 
MAINTAINER Dan Lanier <lanier4@illinois.edu>

ENV SRC_LOC /home

# Clone from github
RUN git clone https://username:password@github.com/KnowEnG-Research/Samples_Clustering_Pipeline.git ${SRC_LOC} 

# Set up working directory
WORKDIR ${SRC_LOC}
 
