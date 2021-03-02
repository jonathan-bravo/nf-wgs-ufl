build: dist/docker
	cp -r ../../modules/ ./modules/
	cp -r ../../bin/ ./bin/
	cp -r ../../main.nf .
	cp -r ../../ufl-germline.nf .
	cp -r ../../ufl-multiqc.nf .
	cp -r ../../ufl-apply_panels.nf .
	cp -r ../../nextflow.config .
	docker build -t nextflow-ufl-germline:1.0.0 .
	rm -rf modules/
	rm -rf bin/
	rm -rf dist/
	rm main.nf
	rm ufl-germline.nf
	rm ufl-multiqc.nf
	rm ufl-apply_panels.nf
	rm nextflow.config
dist/docker:
	curl -fsSLO https://get.docker.com/builds/Linux/x86_64/docker-17.03.1-ce.tgz
	mkdir dist
	tar --strip-components=1 -xvzf docker-17.03.1-ce.tgz -C dist
	rm docker-17.03.1-ce.tgz