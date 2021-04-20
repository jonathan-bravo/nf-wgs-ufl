version ?= 1.0.0

make:
	cp -r ../../Scripts/vcf_report/ClassifyCNV/Resources ./Resources/
	cp -r ../../Scripts/vcf_report/ClassifyCNV/ClassifyCNV.py .
	cp -r ../../Scripts/vcf_report/ClassifyCNV/parse_clingen_tsv.py .
	cp -r ../../Scripts/vcf_report/ClassifyCNV/resources.py .
	cp -r ../../Scripts/vcf_report/varClass.py .
	cp -r ../../Scripts/CNV_json_plot.R .
	cp -r ../../Scripts/g_ranges.py .
	cp -r ../../Reference_Data/PMC-ids.csv.gz .
	aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin 178175818611.dkr.ecr.us-east-1.amazonaws.com
	docker build -t var_class-ufl-germline:${version} .
	docker tag var_class-ufl-germline:${version} 178175818611.dkr.ecr.us-east-1.amazonaws.com/var_class-ufl-germline:${version}
	docker push 178175818611.dkr.ecr.us-east-1.amazonaws.com/var_class-ufl-germline:${version}
	rm -rf Resources/
	rm *.py
	rm *.R
	rm *.gz
