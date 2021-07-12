include ../../common.mk

.PHONY: clean
clean:
	rm -f common/web/templates/index.html
	rm -rf common/web/static
	rm -f common/web/csp-hashes.json

.PHONY: unit-test
unit-test: create-test-db
	PYTHONWARNINGS=ignore:ResourceWarning coverage run \
		--source=app,auth,cli,common,compute,converters,data_anndata,data_common,data_cxg,eb \
		--omit=.coverage,venv \
		-m unittest discover \
		--start-directory ../test/test_czi_hosted/unit/common/apis \
		--top-level-directory ../.. \
		--verbose; test_result=$$?; \
	$(MAKE) clean-test-db; \
	exit $$test_result \


.PHONY: test-db
test-db: create-test-db
	PYTHONWARNINGS=ignore:ResourceWarning coverage run \
		--source=db \
		--omit=.coverage,venv \
		-m unittest discover \
		--start-directory ../test/test_czi_hosted/test_database \
		--top-level-directory ../.. \
		--verbose; test_result=$$?; \
	$(MAKE) clean-test-db; \
	exit $$test_result

.PHONY: create-test-db
create-test-db:
	-docker run -d -p 5432:5432 --name test_db -e POSTGRES_PASSWORD=test_pw postgres

.PHONY: clean-test-db
clean-test-db:
	-docker stop test_db
	-docker rm test_db

.PHONY: test-annotations-performance
test-annotations-performance:
	python ../test/test_czi_hosted/performance/performance_test_annotations_backend.py

.PHONY: test-annotations-scale
test-annotations-scale:
	locust -f ../test/test_czi_hosted/performance/scale_test_annotations.py --headless -u 30 -r 10 --host https://api.cellxgene.dev.single-cell.czi.technology/cellxgene/e/ --run-time 5m 2>&1 | tee locust_dev_stats.txt
