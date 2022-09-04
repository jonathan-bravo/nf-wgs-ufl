version ?= 1.0.0
platform ?= linux/amd64
use_aws ?= true

ifeq (${use_aws}, true)
aws_id ?= $(shell aws sts get-caller-identity --query "Account" --output text)
build: dist/docker
	cp -r ../../modules/ ./modules/
	cp -r ../../bin/ ./bin/
	cp -r ../../workflows/ ./workflows/
	cp -r ../../subworkflows/ ./subworkflows/
	cp -r ../../conf/ ./conf/
	cp -r ../../main.nf .
	cp -r ../../nextflow.config .
	aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin ${aws_id}.dkr.ecr.us-east-1.amazonaws.com
	docker build --platform ${platform} -t nextflow:${version} .
	docker tag nextflow:${version} ${aws_id}.dkr.ecr.us-east-1.amazonaws.com/nextflow-ufl-germline:${version}
	docker push ${aws_id}.dkr.ecr.us-east-1.amazonaws.com/nextflow-ufl-germline:${version}
	rm -rf modules/
	rm -rf bin/
	rm -rf dist/
	rm -rf workflows/
	rm -rf subworkflows/
	rm -rf conf/
	rm *.nf
	rm nextflow.config
dist/docker:
	curl -fsSLO https://get.docker.com/builds/Linux/x86_64/docker-17.03.1-ce.tgz
	mkdir dist
	tar --strip-components=1 -xvzf docker-17.03.1-ce.tgz -C dist
	rm docker-17.03.1-ce.tgz
else
build: dist/docker
	cp -r ../../modules/ ./modules/
	cp -r ../../bin/ ./bin/
	cp -r ../../workflows/ ./workflows/
	cp -r ../../subworkflows/ ./subworkflows/
	cp -r ../../conf/ ./conf/
	cp -r ../../main.nf .
	cp -r ../../nextflow.config .
	docker build --platform ${platform} -t nextflow:${version} .
	rm -rf modules/
	rm -rf bin/
	rm -rf dist/
	rm -rf workflows/
	rm -rf subworkflows/
	rm -rf conf/
	rm *.nf
	rm nextflow.config
dist/docker:
	curl -fsSLO https://get.docker.com/builds/Linux/x86_64/docker-17.03.1-ce.tgz
	mkdir dist
	tar --strip-components=1 -xvzf docker-17.03.1-ce.tgz -C dist
	rm docker-17.03.1-ce.tgz
endif