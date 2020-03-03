include ../common.mk

DATASET := $(if $(DATASET),$(DATASET),../example-dataset/pbmc3k.h5ad)
ANNOTATIONS := $(if $(ANNOTATIONS),$(ANNOTATIONS),../server/test/test_datasets/pbmc3k-annotations.csv)
ANNOTATIONS_FILENAME := $(shell basename $(ANNOTATIONS))

# Packaging
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


# Development convenience methods

.PHONY: server-requirements
server-requirements:
	virtualenv -p python3 ../venv
	source ../venv/bin/activate && \
		pip install -r ../server/requirements-dev.txt && \
		yes | pip uninstall cellxgene || true && \
		pip install -e ..

.PHONY: start-frontend
start-frontend:
	node server/development.js

.PHONY: start-server
start-server:
	cellxgene launch -p $(CXG_SERVER_PORT) $(CXG_OPTIONS) $(DATASET)

.PHONY: backend-dev
backend-dev: server-requirements
	source ../venv/bin/activate && $(MAKE) start-server

.PHONY: backend-dev-anno-ontology
backend-dev-anno-ontology: server-requirements
	CXG_OPTIONS='--annotations-ontology' \
		$(MAKE) backend-dev

.PHONY: test
test:
	node node_modules/jest/bin/jest.js

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
	start_server_and_test '$(MAKE) start-server' $(CXG_SERVER_PORT) '$(MAKE) e2e'


.PHONY: smoke-test-annotations
smoke-test-annotations:
	$(eval TMP_DIR := $(shell  mktemp -d /tmp/cellxgene_XXXXXX))
	cp $(ANNOTATIONS) $(TMP_DIR)/ && \
	start_server_and_test \
		'CXG_OPTIONS="--annotations-input-file $(TMP_DIR)/$(ANNOTATIONS_FILENAME)" $(MAKE) start-server' \
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
