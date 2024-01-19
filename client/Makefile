include ../common.mk

ANNOTATIONS := $(if $(ANNOTATIONS),$(ANNOTATIONS),../test/fixtures/pbmc3k-annotations.csv)
GENE_SETS := $(if $(GENE_SETS),$(GENE_SETS),../test/fixtures/pbmc3k-genesets.csv)
ANNOTATIONS_FILENAME := $(shell basename $(ANNOTATIONS))
GENE_SETS_FILENAME := $(shell basename $(GENE_SETS))

CXG_CONFIG := $(if $(CXG_CONFIG),$(CXG_CONFIG),./__tests__/e2e/test_config.yaml)



# Packaging
.PHONY: clean
clean:
	rm -rf node_modules
	rm -f __tests__/screenshots/*.png

.PHONY: ci
ci:
	npm ci --legacy-peer-deps

.PHONY: install
install:
	npm install

.PHONY: build
WEBPACK_CONFIG ?= configuration/webpack/webpack.config.prod.js
build:
	npm run build $(WEBPACK_CONFIG)

# Development convenience methods
.PHONY: start-frontend
start-frontend:
	node server/development.js

# start an instance of cellxgene and run the end-to-end tests
.PHONY: smoke-test
smoke-test:
	start_server_and_test \
		'CXG_OPTIONS="--config-file $(CXG_CONFIG)" $(MAKE) start-server' \
		$(CXG_SERVER_PORT) \
		'CXG_URL_BASE="http://localhost:$(CXG_SERVER_PORT)" npm run e2e -- --verbose false'

# start an instance of cellxgene and run the end-to-end annotations tests
.PHONY: smoke-test-annotations
smoke-test-annotations:
	$(eval TMP_DIR := $(shell  mktemp -d /tmp/cellxgene_XXXXXX))
	cp $(ANNOTATIONS) $(TMP_DIR)/ && \
	cp $(GENE_SETS) $(TMP_DIR)/ && \
	start_server_and_test \
		'CXG_OPTIONS="--annotations-file $(TMP_DIR)/$(ANNOTATIONS_FILENAME) --gene-sets-file $(TMP_DIR)/$(GENE_SETS_FILENAME)" $(MAKE) start-server' \
		$(CXG_SERVER_PORT) \
		'CXG_URL_BASE="http://localhost:$(CXG_SERVER_PORT)" npm run e2e-annotations -- --verbose false'
	rm -rf $(TMP_DIR)

.PHONY: unit-test
unit-test:
	node node_modules/jest/bin/jest.js --testPathIgnorePatterns e2e

# pass remaining commands through to npm run
%:
	npm run $(*)
