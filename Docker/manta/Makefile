version ?= 1.0.0
platform ?= linux/amd64
use_aws ?= true

ifeq (${use_aws}, true)
aws_id ?= $(shell aws sts get-caller-identity --query "Account" --output text)
make:
	aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin ${aws_id}.dkr.ecr.us-east-1.amazonaws.com
	docker build --platform ${platform} -t manta:${version} .
	docker tag manta:${version} ${aws_id}.dkr.ecr.us-east-1.amazonaws.com/manta-ufl-germline:${version}
	docker push ${aws_id}.dkr.ecr.us-east-1.amazonaws.com/manta-ufl-germline:${version}
else
make:
	docker build --platform ${platform} -t manta:${version} .
endif