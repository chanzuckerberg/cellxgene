include common.mk

BUILDDIR := build
CLIENTBUILD := $(BUILDDIR)/client
CZIHOSTEDBUILD := $(BUILDDIR)/backend/czi_hosted
SERVERBUILD := $(BUILDDIR)/backend/server
CLEANFILES :=  $(BUILDDIR)/ client/build build dist cellxgene.egg-info

PART ?= patch

# CLEANING
.PHONY: clean
clean: clean-lite clean-czi-hosted clean-server clean-client

# cleaning the client's node_modules is the longest one, so we avoid that if possible
.PHONY: clean-lite
clean-lite:
	rm -rf $(CLEANFILES)

.PHONY: clean-client
clean-client:
	cd client && $(MAKE) clean

.PHONY: clean-server
clean-server:
	cd backend/server && $(MAKE) clean

.PHONY: clean-czi-hosted
clean-czi-hosted:
	cd backend/czi_hosted && $(MAKE) clean

# BUILDING PACKAGE

.PHONY: build-client
build-client:
	cd client && $(MAKE) ci build

.PHONY: build
build: clean build-client
	git ls-files backend/server/ | grep -v 'backend/server/test/' | cpio -pdm $(BUILDDIR)
	cp -r client/build/  $(CLIENTBUILD)
	$(call copy_client_assets,$(CLIENTBUILD),$(SERVERBUILD))
	cp backend/__init__.py $(BUILDDIR)
	cp backend/__init__.py $(BUILDDIR)/backend
	cp -r backend/common $(BUILDDIR)/backend/common
	cp MANIFEST.in README.md setup.cfg setup.py $(BUILDDIR)

.PHONY: build-czi-hosted
build-czi-hosted: clean build-client
	git ls-files backend/czi_hosted/ | grep -v 'backend/czi_hosted/test/' | cpio -pdm $(BUILDDIR)
	cp -r client/build/  $(CLIENTBUILD)
	$(call copy_client_assets,$(CLIENTBUILD),$(CZIHOSTEDBUILD))
	cp -r backend/common $(BUILDDIR)/backend/common
	cp backend/__init__.py $(BUILDDIR)
	cp backend/__init__.py $(BUILDDIR)/backend
	cp MANIFEST_hosted.in README.md setup.cfg setup_hosted.py $(BUILDDIR)
	mv $(BUILDDIR)/setup_hosted.py $(BUILDDIR)/setup.py
	mv $(BUILDDIR)/MANIFEST_hosted.in $(BUILDDIR)/MANIFEST.in

# If you are actively developing in the server folder use this, dirties the source tree
.PHONY: build-for-server-dev
build-for-server-dev: clean-server build-client
	$(call copy_client_assets,client/build,backend/server)

.PHONY: build-for-czi-hosted-dev
build-for-czi-hosted-dev: clean-czi-hosted build-client
	$(call copy_client_assets,client/build,backend/czi_hosted)

.PHONY: copy-client-assets
copy-client-assets:
	$(call copy_client_assets,client/build,backend/server)

.PHONY: copy-client-assets-czi-hosted
copy-client-assets-czi-hosted:
	$(call copy_client_assets,client/build,backend/czi_hosted)

# TESTING
.PHONY: test
test: unit-test smoke-test

.PHONY: unit-test
unit-test: unit-test-server unit-test-client unit-test-common

.PHONY: test-server
test-server: unit-test-server smoke-test

.PHONY: test-czi-hosted
test-czi-hosted: unit-test-czi-hosted smoke-test

.PHONY: unit-test-client
unit-test-client:
	cd client && $(MAKE) unit-test

.PHONY: unit-test-czi-hosted
unit-test-czi-hosted:
	cd backend/czi_hosted && $(MAKE) unit-test

.PHONY: unit-test-server
unit-test-server:
	cd backend/server && $(MAKE) unit-test

.PHONY: unit-test-common
unit-test-common:
	cd backend/common && $(MAKE) unit-test

.PHONY: smoke-test
smoke-test:
	cd client && $(MAKE) smoke-test

.PHONY: smoke-test-annotations
smoke-test-annotations:
	cd client && $(MAKE) smoke-test-annotations

.PHONY: test-db
test-db:
	cd backend/czi_hosted && $(MAKE) test-db

# FORMATTING CODE

.PHONY: fmt
fmt: fmt-client fmt-py

.PHONY: fmt-client
fmt-client:
	cd client && $(MAKE) fmt

.PHONY: fmt
fmt-py:
	black .

.PHONY: lint
lint: lint-servers lint-client

.PHONY: lint-servers
lint-servers: lint-server lint-czi-hosted-server

.PHONY: lint-server
lint-server: fmt-py
	flake8 backend/server --per-file-ignores='backend/test/fixtures/dataset_config_outline.py:F821 backend/test/fixtures/server_config_outline.py:F821 backend/server/test/performance/scale_test_annotations.py:E501'

.PHONY: lint-czi-hosted-server
lint-czi-hosted-server: fmt-py
	flake8 backend/czi_hosted --per-file-ignores='backend/test/fixtures/czi_hosted_dataset_config_outline.py:F821 backend/test/fixtures/czi_hosted_server_config_outline.py:F821 backend/test/performance/scale_test_annotations.py:E501'

.PHONY: lint-client
lint-client:
	cd client && $(MAKE) lint

# CREATING DISTRIBUTION RELEASE

.PHONY: pydist
pydist: build
	cd $(BUILDDIR); python setup.py sdist -d ../dist
	@echo "done"

.PHONY: pydist
pydist-czi-hosted: build-czi-hosted
	cd $(BUILDDIR); python setup.py sdist -d ../dist
	@echo "done"


# RELEASE HELPERS

# Create new version to commit to main
.PHONY: create-release-candidate
create-release-candidate: dev-env bump-version clean-lite gen-package-lock
	@echo "Version bumped part:$(PART) and client built. Ready to commit and push"

# Bump the release candidate version if needed (i.e. the previous release candidate had errors).
.PHONY: recreate-release-candidate
recreate-release-candidate: dev-env bump-release-candidate clean-lite gen-package-lock
	@echo "Version bumped part:$(PART) and client built. Ready to commit and push"

# Build dist and release to Test PyPI
.PHONY: release-candidate-to-test-pypi
release-candidate-to-test-pypi: dev-env pydist twine
	@echo "Dist built and uploaded to test.pypi.org"
	@echo "Test the install:"
	@echo "    make install-release-test"

# Build final dist (gets rid of the rc tag) and release final candidate to TestPyPI
.PHONY: release-final-to-test-pypi
release-final-to-test-pypi: dev-env bump-release clean-lite gen-package-lock pydist twine
	@echo "Final release dist built and uploaded to test.pypi.org"
	@echo "Test the install:"
	@echo "    make install-release-test"

.PHONY: release-final
release-final: twine-prod
	@echo "Release uploaded to pypi.org"

# DANGER: releases directly to prod
# use this if you accidently burned a test release version number,
.PHONY: release-directly-to-prod
release-directly-to-prod: dev-env pydist twine-prod
	@echo "Dist built and uploaded to pypi.org"
	@echo "Test the install:"
	@echo "    make install-release"

.PHONY: dev-env
dev-env: dev-env-client dev-env-server

.PHONY: dev-env-client
dev-env-client:
	cd client && $(MAKE) ci

.PHONY: dev-env-server
dev-env-server:
	pip install -r backend/server/requirements-dev.txt

.PHONY: dev-env-czi-hosted
dev-env-czi-hosted:
	pip install -r backend/czi_hosted/requirements-dev.txt
# Set PART=[major, minor, patch] as param to make bump.
# This will create a release candidate. (i.e. 0.16.1 -> 0.16.2-rc.0 for a patch bump)
.PHONY: bump-version
bump-version:
	bumpversion --config-file .bumpversion.cfg $(PART)

# Increments the release candidate version (i.e. 0.16.2-rc.1 -> 0.16.2-rc.2)
.PHONY: bump-release-candidate
bump-release-candidate:
	bumpversion --config-file .bumpversion.cfg prerelversion --allow-dirty

# Finalizes the release candidate by removing the release candidate tag (i.e. 0.16.2-rc.2 -> 0.16.2).
.PHONY: bump-release
bump-release:
	bumpversion --config-file .bumpversion.cfg prerel --allow-dirty

.PHONY: twine
twine:
	twine upload --repository-url https://test.pypi.org/legacy/ dist/*

.PHONY: twine-prod
twine-prod:
	twine upload dist/*

# quicker than re-building client
.PHONY: gen-package-lock
gen-package-lock:
	cd client && $(MAKE) install


# INSTALL

# install from source tree for development
.PHONY: install-dev
install-dev: uninstall
	pip install -e .

# install from test.pypi to test your release
.PHONY: install-release-test
install-release-test: uninstall
	pip install --no-cache-dir --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple cellxgene
	@echo "Installed cellxgene from test.pypi.org, now run and smoke test"

# install from pypi to test your release
.PHONY: install-release
install-release: uninstall
	pip install --no-cache-dir cellxgene
	@echo "Installed cellxgene from pypi.org"

# install from dist
.PHONY: install-dist
install-dist: uninstall
	pip install dist/cellxgene*.tar.gz

.PHONY: uninstall
uninstall:
	pip uninstall -y cellxgene || :

