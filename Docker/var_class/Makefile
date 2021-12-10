version ?= 1.0.0
platform ?= linux/amd64

make:
	cp -r ../../Scripts/vcf_report/ClassifyCNV/Resources ./Resources/
	cp -r ../../Scripts/vcf_report/ClassifyCNV/ClassifyCNV.py .
	cp -r ../../Scripts/vcf_report/ClassifyCNV/parse_clingen_tsv.py .
	cp -r ../../Scripts/vcf_report/ClassifyCNV/resources.py .
	cp -r ../../Scripts/vcf_report/reporting.py .
	cp -r ../../Scripts/json_to_csv.py .
	aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin 178175818611.dkr.ecr.us-east-1.amazonaws.com
	docker buildx build --platform ${platform} -t reporting:${version} .
	docker tag reporting:${version} 178175818611.dkr.ecr.us-east-1.amazonaws.com/var_class-ufl-germline:${version}
	docker push 178175818611.dkr.ecr.us-east-1.amazonaws.com/var_class-ufl-germline:${version}
	rm -rf Resources/
	rm *.py
