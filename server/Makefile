include ../common.mk

.PHONY: clean
clean:
	rm -f common/web/templates/index.html
	rm -rf common/web/static
	rm -f common/web/csp-hashes.json


.PHONY: test-annotations-performance
test-annotations-performance:
	python ../test/performance/performance_test_annotations.py
