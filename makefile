BUILDDIR := build
CLIENTBUILD := $(BUILDDIR)/client
SERVERBUILD := $(BUILDDIR)/server
CLEANFILES :=  $(BUILDDIR)/ client/build dist cellxgene.egg-info


install:
	pip install -e .

uninstall:
	pip uninstall cellxgene

build: clean-lite build-server
	@echo "done"

pydist: clean build-server
	python setup.py sdist
	@echo "done"

build-server: build-client
	mkdir -p $(SERVERBUILD)
	cp -r server/* $(SERVERBUILD)
	cp -r client/build/  $(CLIENTBUILD)
	mkdir -p $(SERVERBUILD)/app/web/static/img
	cp $(CLIENTBUILD)/index.html $(SERVERBUILD)/app/web/templates/
	cp -r $(CLIENTBUILD)/static $(SERVERBUILD)/app/web/
	cp $(CLIENTBUILD)/favicon.png $(SERVERBUILD)/app/web/static/img
	cp $(CLIENTBUILD)/service-worker.js $(SERVERBUILD)/app/web/static/js/

build-client: 
	npm install --prefix client/ client
	npm run  --prefix client build

clean : clean-lite
	rm -rf client/node_modules

clean-lite :
	rm -rf $(CLEANFILES)

clean-server:
	rm -f server/app/web/templates/index.html
	rm -rf server/app/web/static

build-for-server-dev: clean-server
	mkdir -p server/app/web/static/img
	cp client/build/index.html server/app/web/templates/
	cp -r client/build/static server/app/web/
	cp client/build/favicon.png server/app/web/static/img
	cp client/build/service-worker.js server/app/web/static/js/

charlotte:
	rm -rf dist cellxgene.egg-info
	yes | pip uninstall cellxgene || true
	pip install -e .
	cellxgene launch example-dataset/pbmc3k.h5ad


.PHONY: pydist build build-server build-client clean clean-lite clean-server build-for-server-dev