version ?= 1.0.0
platform ?= linux/amd64

build: dist/docker
	cp -r ../../modules/ ./modules/
	cp -r ../../bin/ ./bin/
	cp -r ../../workflows/ ./workflows/
	cp -r ../../*.nf .
	cp -r ../../nextflow.config .
	docker buildx build --platform ${platform} -t local_nextflow:${version} .
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
