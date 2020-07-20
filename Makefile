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

clean-%:
	cd $(*) && $(MAKE) clean


# BUILDING PACKAGE

.PHONY: build-client
build-client:
	cd client && $(MAKE) ci build

.PHONY: build
build: clean build-client
	git ls-files server/ | grep -v 'server/test/' | cpio -pdm $(BUILDDIR)
	cp -r client/build/  $(CLIENTBUILD)
	$(call copy_client_assets,$(CLIENTBUILD),$(SERVERBUILD))
	cp MANIFEST.in README.md setup.cfg setup.py $(BUILDDIR)

# If you are actively developing in the server folder use this, dirties the source tree
.PHONY: build-for-server-dev
build-for-server-dev: clean-server build-client
	$(call copy_client_assets,client/build,server)

.PHONY: copy-client-assets
copy-client-assets:
	$(call copy_client_assets,client/build,server)

# TESTING
.PHONY: test
test: unit-test smoke-test

.PHONY: unit-test
unit-test: unit-test-server unit-test-client

unit-test-%:
	cd $(*) && $(MAKE) unit-test

.PHONY: smoke-test
smoke-test:
	cd client && $(MAKE) smoke-test

.PHONY: smoke-test-annotations
smoke-test-annotations:
	cd client && $(MAKE) smoke-test-annotations

# FORMATTING CODE

.PHOHY: fmt
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
lint-server:
	flake8 server

.PHONY: lint-client
lint-client:
	cd client && $(MAKE) lint

# CREATING DISTRIBUTION RELEASE

.PHONY: pydist
pydist: build
	cd $(BUILDDIR); python setup.py sdist -d ../dist
	@echo "done"


# RELEASE HELPERS

# create new version to commit to main
.PHONY: release-stage-1
release-stage-1: dev-env bump clean-lite gen-package-lock
	@echo "Version bumped part:$(PART) and client built. Ready to commit and push"

# build dist and release to dev pypi
.PHONY: release-stage-2
release-stage-2: dev-env pydist twine
	@echo "Dist built and uploaded to test.pypi.org"
	@echo "Test the install:"
	@echo "    make install-release-test"
	@echo "Then upload to Pypi prod:"
	@echo "    make twine-prod"

.PHONY: release-stage-final
release-stage-final: twine-prod
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

# give PART=[major, minor, part] as param to make bump
.PHONY: bump
bump:
	bumpversion --config-file .bumpversion.cfg $(PART)

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

