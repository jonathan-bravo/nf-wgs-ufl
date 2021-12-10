make:
	aws s3 sync ../pipeline_webpage/ s3://hakmonkey-ufl-germline-landing-page/
	aws cloudfront create-invalidation --distribution-id E24KGUSZG5V2YM --invalidation-batch file://landing-page_refresh.json