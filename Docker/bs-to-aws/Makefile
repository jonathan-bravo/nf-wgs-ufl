version ?= 1.0
platform ?= linux/amd64

make:
	cp ../../Scripts/bs-to-aws.sh .
	cp -r ~/.basespace .
	aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin 178175818611.dkr.ecr.us-east-1.amazonaws.com
	docker build --platform ${platform} -t bs-to-aws:${version} .
	docker tag bs-to-aws:${version} 178175818611.dkr.ecr.us-east-1.amazonaws.com/bs-to-aws:${version}
	docker push 178175818611.dkr.ecr.us-east-1.amazonaws.com/bs-to-aws:${version}
	rm bs-to-aws.sh
	rm -rf .basespace/
