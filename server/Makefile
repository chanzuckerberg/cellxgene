include ../common.mk

.PHONY: clean
clean:
	rm -f common/web/templates/index.html
	rm -rf common/web/static

.PHONY: unit-test
unit-test:
	pytest -s test
