include common.mk

BUILDDIR := build
CLIENTBUILD := $(BUILDDIR)/client
SERVERBUILD := $(BUILDDIR)/server
CLEANFILES :=  $(BUILDDIR)/ client/build build dist cellxgene.egg-info

PART ?= patch

# CLEANING
.PHONY: clean
clean: clean-lite clean-server clean-client

# cleaning the client's node_modules is the longest one, so we avoid that if possible
.PHONY: clean-lite
clean-lite:
	rm -rf $(CLEANFILES)

.PHONY: clean-client
clean-client:
	cd client && $(MAKE) clean

.PHONY: clean-server
clean-server:
	cd server && $(MAKE) clean


# BUILDING PACKAGE

.PHONY: build-client
build-client:
	cd client && $(MAKE) ci build

.PHONY: build
build: clean build-client
	git ls-files server/ | cpio -pdm $(BUILDDIR)
	cp -r client/build/  $(CLIENTBUILD)
	$(call copy_client_assets,$(CLIENTBUILD),$(SERVERBUILD))
	cp MANIFEST.in README.md setup.cfg setup.py $(BUILDDIR)

# If you are actively developing in the server folder use this, dirties the source tree
.PHONY: build-for-server-dev
build-for-server-dev: clean-server build-client copy-client-assets

.PHONY: copy-client-assets
copy-client-assets:
	$(call copy_client_assets,client/build,server)


# TESTING
.PHONY: test
test: unit-test smoke-test

.PHONY: unit-test
unit-test: unit-test-server unit-test-client

.PHONY: test-server
test-server: unit-test-server smoke-test

.PHONY: unit-test-client
unit-test-client:
	cd client && $(MAKE) unit-test

.PHONY: unit-test-server
unit-test-server:
	PYTHONWARNINGS=ignore:ResourceWarning coverage run \
		--source=server \
		--omit=.coverage,venv \
		-m unittest discover \
		--start-directory test/unit \
		--verbose; test_result=$$?; \
	exit $$test_result \

.PHONY: smoke-test
smoke-test:
	cd client && $(MAKE) smoke-test

.PHONY: smoke-test-annotations
smoke-test-annotations:
	cd client && $(MAKE) smoke-test-annotations

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
lint: lint-server lint-client


.PHONY: lint-server
lint-server: fmt-py
	flake8 server --per-file-ignores='test/fixtures/dataset_config_outline.py:F821 test/fixtures/server_config_outline.py:F821 test/performance/scale_test_annotations.py:E501'

.PHONY: lint-client
lint-client:
	cd client && $(MAKE) lint

# CREATING DISTRIBUTION RELEASE

.PHONY: pydist
pydist: build
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
	pip install -r server/requirements-dev.txt

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

