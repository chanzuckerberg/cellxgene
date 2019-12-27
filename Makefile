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

.PHONY: build
build: clean build-server
	@echo "done"

.PHONY: build-client
build-client:
	cd client && $(MAKE) install build

.PHONY: build-server
build-server: build-client
	mkdir -p $(SERVERBUILD)
	cp -r server/* $(SERVERBUILD)
	cp -r client/build/  $(CLIENTBUILD)
	mkdir -p $(SERVERBUILD)/app/web/static/img
	mkdir -p $(SERVERBUILD)/app/web/templates/
	cp $(CLIENTBUILD)/index.html $(SERVERBUILD)/app/web/templates/
	cp -r $(CLIENTBUILD)/static $(SERVERBUILD)/app/web/
	cp $(CLIENTBUILD)/favicon.png $(SERVERBUILD)/app/web/static/img
	cp $(CLIENTBUILD)/service-worker.js $(SERVERBUILD)/app/web/static/js/
	cp MANIFEST.in README.md setup.cfg setup.py $(BUILDDIR)

# If you are actively developing in the server folder use this, dirties the source tree
.PHONY: build-for-server-dev
build-for-server-dev: clean-server build-client
	mkdir -p server/app/web/static/img
	mkdir -p server/app/web/static/js
	mkdir -p server/app/web/templates/
	cp client/build/index.html server/app/web/templates/
	cp -r client/build/static server/app/web/
	cp client/build/favicon.png server/app/web/static/img
	cp client/build/service-worker.js server/app/web/static/js/


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

# FORMATTING CODE

.PHOHY: fmt
fmt: fmt-client fmt-py

fmt-client:
	cd client && $(MAKE) fmt

fmt-py:
	black .

.PHONY: lint
lint:
	flake8 server


# CREATING DISTRIBUTION RELEASE

.PHONY: pydist
pydist: build
	cd $(BUILDDIR); python setup.py sdist -d ../dist
	@echo "done"


# RELEASE HELPERS

# create new version to commit to master
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
dev-env:
	cd client && $(MAKE) install
	pip install -r server/requirements-dev.txt

.PHONY: gui-env
gui-env: dev-env
	pip install -r server/requirements-gui.txt

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

# setup.py sucks when you have your library in a separate folder, adding these in to help setup envs

# install from build directory
.PHONY: install
install: uninstall
	cd $(BUILDDIR); pip install -e .

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


# GUI

.PHONY: build-assets
build-assets:
	pyside2-rcc server/gui/cellxgene.qrc -o server/gui/cellxgene_rc.py

.PHONY: gui-spec-osx
gui-spec-osx: clean-lite gui-env
	pip install -e .[gui]
	pyi-makespec -D -w --additional-hooks-dir server/gui/ -n cellxgene  --add-binary='/System/Library/Frameworks/Tk.framework/Tk':'tk' --add-binary='/System/Library/Frameworks/Tcl.framework/Tcl':'tcl'  --add-data server/app/web/templates/:server/app/web/templates/ --add-data server/app/web/static/:server/app/web/static/ --icon server/gui/images/cxg_icons.icns server/gui/main.py
	mv cellxgene.spec cellxgene-osx.spec

.PHONY: gui-spec-windows
gui-spec-windows: clean-lite dev-env
	pip install -e .[gui]
	pyi-makespec -D -w --additional-hooks-dir server/gui/ -n cellxgene --add-data server/app/web/templates;server/app/web/templates --add-data server/app/web/static;server/app/web/static --icon server/gui/images/icon.ico server/gui/main.py
	mv cellxgene.spec cellxgene-windows.spec

.PHONY: gui-build-osx
gui-build-osx: clean-lite
	pyinstaller --clean cellxgene-osx.spec

.PHONY: gui-build-windows
gui-build-windows: clean-lite
	pyinstaller --clean cellxgene-windows.spec

