PROJECT_ROOT := $(shell git rev-parse --show-toplevel)
PATH := $(PATH):$(PROJECT_ROOT)/scripts

SHELL := env PATH='$(PATH)' /bin/bash

ifeq ($(shell which jq),)
$(error Please install jq using "apt-get install jq" or "brew install jq")
endif

define GetValueFromJson
$(shell jq -r '.$(1)' "$(PROJECT_ROOT)/environment.default.json")
endef

# https://stackoverflow.com/a/14777895/9587410
ifeq ($(shell uname),Darwin)     # is Windows_NT on XP, 2000, 7, Vista, 10...
    IS_DARWIN := "true"
endif

export CELLXGENE_COMMIT := $(shell git rev-parse --short HEAD)
export CXG_SERVER_PORT := $(call GetValueFromJson,CXG_SERVER_PORT)
export CXG_CLIENT_PORT := $(call GetValueFromJson,CXG_CLIENT_PORT)
export JEST_ENV := $(call GetValueFromJson,JEST_ENV)
export DATASET := $(call GetValueFromJson,DATASET)
export CXG_OPTIONS := $(call GetValueFromJson,CXG_OPTIONS)

.PHONY: start-server
start-server:
	cellxgene launch -p $(CXG_SERVER_PORT) $(CXG_OPTIONS) $(DATASET)

# copy the client assests to a location known to the server
# $(1) is the source of the client assets
# $(2) is the destination
define copy_client_assets
	mkdir -p $(2)/common/web/static/assets
	mkdir -p $(2)/common/web/templates/
	cp $(1)/index.html $(2)/common/web/templates/
	cp -r $(1)/static $(2)/common/web/
	cp $(1)/csp-hashes.json $(2)/common/web/
endef

