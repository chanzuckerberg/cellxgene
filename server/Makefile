.PHONY: fmt
fmt:
	yapf -ir .

.PHONY: lint
lint:
	flake8

.PHONY: clean
clean:
	rm -f app/web/templates/index.html
	rm -rf app/web/static

.PHONY: unit-test
unit-test:
	pytest -s test
