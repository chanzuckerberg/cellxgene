include ../../common.mk

.PHONY: clean
clean:
	rm -f artifact.zip
	rm -rf artifact.dir


# Build the ElasticBeanstalk configuration and deployment bundle,
# such that deployment can be done with a simple `eb deploy`.
# Presumes that a top-level `make build-client` has been done to 
# create the client static assets.

.PHONY: build
build: clean
	mkdir artifact.dir; \
	(cd ../.. ; \
	    git ls-files server/ | cpio -pdm server/eb/artifact.dir ; ); \
	$(call copy_client_assets,../../client/build,artifact.dir/server) ; \
	cp app.py artifact.dir/application.py; \
	cp ../requirements.txt artifact.dir; \
	cp -r .ebextensions artifact.dir; \
	(cd artifact.dir; zip -r ../artifact.zip . --exclude server/test/\* server/eb/\* ); \
	if ! grep artifact.zip .elasticbeanstalk/config.yml ; then \
	    mkdir -p .elasticbeanstalk ; cat config_deploy.yaml >> .elasticbeanstalk/config.yml ; fi

