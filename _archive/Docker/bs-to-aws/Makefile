version ?= 1.0
platform ?= linux/amd64
use_aws ?= true

ifeq (${use_aws}, true)
aws_id ?= $(shell aws sts get-caller-identity --query "Account" --output text)
make:
	cp ../../Scripts/bs-to-aws.sh .
	cp -r ~/.basespace .
	aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin ${aws_id}.dkr.ecr.us-east-1.amazonaws.com
	docker build --platform ${platform} -t bs-to-aws:${version} .
	docker tag bs-to-aws:${version} ${aws_id}.dkr.ecr.us-east-1.amazonaws.com/bs-to-aws:${version}
	docker push ${aws_id}.dkr.ecr.us-east-1.amazonaws.com/bs-to-aws:${version}
	rm bs-to-aws.sh
	rm -rf .basespace/
else
make:
	cp ../../Scripts/bs-to-aws.sh .
	cp -r ~/.basespace .
	docker build --platform ${platform} -t bs-to-aws:${version} .
	rm bs-to-aws.sh
	rm -rf .basespace/
endif