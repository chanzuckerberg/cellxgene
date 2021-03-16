
.PHONY: unit-test
unit-test:
	PYTHONWARNINGS=ignore:ResourceWarning coverage run \
		--omit=.coverage,data_common/fbs/NetEncoding,venv \
		-m unittest discover \
		--start-directory ../test/unit/test_common \
		--top-level-directory ../../ \
		--verbose; test_result=$$?; \
	exit $$test_result \
