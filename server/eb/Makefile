include ../../common.mk

.PHONY: clean
clean:
	rm -f artifact.zip
	rm -rf artifact.dir


# Build the ElasticBeanstalk configuration and deployment bundle,
# such that deployment can be done with a simple `eb deploy`.
# Presumes that a top-level `make build-client` has been done to 
# create the client static assets.

cwd := $(shell pwd)

.PHONY: build
build: clean
	mkdir artifact.dir; \
	(cd ../.. ; \
	  git ls-files server/ | cpio -pdm $(cwd)/artifact.dir ; ); \
	$(call copy_client_assets,../../client/build,artifact.dir/server) ; \
	set -e ; \
	cp app.py artifact.dir/application.py; \
	cp ../requirements.txt artifact.dir; \
	cp -r .ebextensions artifact.dir; \
	if [ -d customize ] ; then \
	  if [ -f customize/config.yaml ] ; then \
	     cp customize/config.yaml artifact.dir;  \
	  fi ; \
	  if [ -f customize/requirements.txt ] ; then \
	     pip install requirements-parser ; \
	     pip install packaging ; \
	     python3 check_requirements.py ../requirements.txt customize/requirements.txt; \
	     cp customize/requirements.txt artifact.dir;  \
	  fi ; \
	  if [ -d customize/deploy ] ; then \
		 mkdir -p artifact.dir/server/common/web/static/cellxgene; \
		 cp -r customize/deploy artifact.dir/server/common/web/static/cellxgene;  \
	  fi; \
	  if [ -d customize/inline_scripts ] ; then \
	      cp -r customize/inline_scripts/* artifact.dir/server/common/web/templates;  \
	  fi; \
	  if [ -d customize/ebextensions ] ; then \
	      cp -r customize/ebextensions/* artifact.dir/.ebextensions;  \
	  fi; \
	fi; \
	if [ -d customize/plugins ] ; then \
	  cp -r customize/plugins artifact.dir/server/;  \
	fi; \
	(cd artifact.dir; \
	  cp -r server/common/web/static static; \
	  zip -r ../artifact.zip . --exclude server/test/\* server/eb/\* ; ); \
	if ! grep artifact.zip .elasticbeanstalk/config.yml ; then \
	    mkdir -p .elasticbeanstalk ; cat config_deploy.yaml >> .elasticbeanstalk/config.yml ; fi

