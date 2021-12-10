version ?= 1.0.0
platform ?= linux/amd64

build: dist/docker
	cp -r ../../modules/ ./modules/
	cp -r ../../bin/ ./bin/
	cp -r ../../workflows/ ./workflows/
	cp -r ../../main.nf .
	cp -r ../../nextflow.config .
	aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin 178175818611.dkr.ecr.us-east-1.amazonaws.com
	docker buildx build --platform ${platform} -t nextflow:${version} .
	docker tag nextflow:${version} 178175818611.dkr.ecr.us-east-1.amazonaws.com/nextflow-ufl-germline:${version}
	docker push 178175818611.dkr.ecr.us-east-1.amazonaws.com/nextflow-ufl-germline:${version}
	rm -rf modules/
	rm -rf bin/
	rm -rf dist/
	rm -rf workflows/
	rm *.nf
	rm nextflow.config
dist/docker:
	curl -fsSLO https://get.docker.com/builds/Linux/x86_64/docker-17.03.1-ce.tgz
	mkdir dist
	tar --strip-components=1 -xvzf docker-17.03.1-ce.tgz -C dist
	rm docker-17.03.1-ce.tgz
