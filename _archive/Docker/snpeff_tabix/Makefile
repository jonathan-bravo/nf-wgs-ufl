version ?= 1.0.0
platform ?= linux/amd64
use_aws ?= true

ifeq (${use_aws}, true)
aws_id ?= $(shell aws sts get-caller-identity --query "Account" --output text)
make:
	aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin 178175818611${aws_id}.dkr.ecr.us-east-1.amazonaws.com
	docker build --platform ${platform} -t snpeff:${version} .
	docker tag snpeff:${version} ${aws_id}.dkr.ecr.us-east-1.amazonaws.com/snpeff_tabix-ufl-germline:${version}
	docker push ${aws_id}.dkr.ecr.us-east-1.amazonaws.com/snpeff_tabix-ufl-germline:${version}
else
make:
	docker build --platform ${platform} -t snpeff:${version} .
endif