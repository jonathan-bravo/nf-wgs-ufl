version ?= 1.0.0
platform ?= linux/amd64

make:
	aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin 178175818611.dkr.ecr.us-east-1.amazonaws.com
	docker buildx build --platform ${platform} -t multiqc:${version} .
	docker tag multiqc:${version} 178175818611.dkr.ecr.us-east-1.amazonaws.com/multiqc-ufl-germline:${version}
	docker push 178175818611.dkr.ecr.us-east-1.amazonaws.com/multiqc-ufl-germline:${version}
