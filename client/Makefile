.PHONY: clean
clean:
	rm -rf node_modules

.PHONY: ci
ci:
	npm ci client

.PHONY: install
install:
	npm install client

.PHONY: build
build:
	npm run build

# pass remaining commands through to npm run
%:
	npm run $(*)

