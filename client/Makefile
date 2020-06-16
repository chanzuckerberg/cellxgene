include ../common.mk

ANNOTATIONS := $(if $(ANNOTATIONS),$(ANNOTATIONS),../server/test/test_datasets/pbmc3k-annotations.csv)
ANNOTATIONS_FILENAME := $(shell basename $(ANNOTATIONS))

# Packaging
.PHONY: clean
clean:
	rm -rf node_modules

.PHONY: ci
ci:
	npm ci

.PHONY: install
install:
	npm install

.PHONY: build
WEBPACK_CONFIG ?= configuration/webpack/webpack.config.prod.js
build:
	npm run build $(WEBPACK_CONFIG)

# Formatting code

.PHONY: lint
lint:
	npx eslint ./src/

# Development convenience methods
.PHONY: start-frontend
start-frontend:
	node server/development.js

.PHONY: e2e
e2e:
	node node_modules/jest/bin/jest.js \
		--verbose false \
		--config __tests__/e2e/e2eJestConfig.json \
		e2e/e2e.test.js

.PHONY: e2e-annotations
e2e-annotations:
	node node_modules/jest/bin/jest.js \
		--verbose false \
		--config __tests__/e2e/e2eJestConfig.json \
		e2e/e2eAnnotations.test.js

.PHONY: smoke-test
smoke-test:
	start_server_and_test \
		'CXG_OPTIONS="--disable-annotations" $(MAKE) start-server' \
		$(CXG_SERVER_PORT) \
		'$(MAKE) e2e'


.PHONY: smoke-test-annotations
smoke-test-annotations:
	$(eval TMP_DIR := $(shell  mktemp -d /tmp/cellxgene_XXXXXX))
	cp $(ANNOTATIONS) $(TMP_DIR)/ && \
	start_server_and_test \
		'CXG_OPTIONS="--annotations-file $(TMP_DIR)/$(ANNOTATIONS_FILENAME)" $(MAKE) start-server' \
		$(CXG_SERVER_PORT) \
		'$(MAKE) e2e-annotations'
	rm -rf $(TMP_DIR)

.PHONY: unit-test
unit-test:
	node node_modules/jest/bin/jest.js \
	    --testPathIgnorePatterns e2e

# pass remaining commands through to npm run
%:
	npm run $(*)
