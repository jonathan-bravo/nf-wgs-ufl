version ?= 1.0.0
platform ?= linux/amd64

make:
	aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin 178175818611.dkr.ecr.us-east-1.amazonaws.com
	docker build --platform ${platform} -t samtools:${version} .
	docker tag samtools:${version} 178175818611.dkr.ecr.us-east-1.amazonaws.com/samtools-ufl-germline:${version}
	docker push 178175818611.dkr.ecr.us-east-1.amazonaws.com/samtools-ufl-germline:${version}
