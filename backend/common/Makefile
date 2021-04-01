
.PHONY: unit-test
unit-test:
	PYTHONWARNINGS=ignore:ResourceWarning coverage run \
		--source=fbs,utils \
		--omit=.coverage,data_common/fbs/NetEncoding,venv \
		-m unittest discover \
		--start-directory ../test/test_common/unit \
		--top-level-directory ../../ \
		--verbose; test_result=$$?; \
	exit $$test_result \
