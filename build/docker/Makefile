CUR_DATE=`date +'%m.%d.%Y'`

.PHONY: build_docker_image

build_docker_image:
	docker build -t knowengdev/samples_clustering_pipeline:$(CUR_DATE) .

push_to_dockerhub:
	docker push knowengdev/samples_clustering_pipeline
