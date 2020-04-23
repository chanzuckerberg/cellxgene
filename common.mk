PROJECT_ROOT := $(shell git rev-parse --show-toplevel)
PATH := $(PATH):$(PROJECT_ROOT)/scripts

SHELL := env PATH='$(PATH)' /bin/bash

# get_or_else_dev_env_default
# - If a variable is defined, return its value
# - Else return the default value from environment.dev
define get_or_else_dev_env_default
$(if $($(1)),$($(1)),$(shell VAR=$$(sed -n 's/$(1)=\(.*\)/\1/p' $(PROJECT_ROOT)/environment.default); eval "echo \"$$VAR\""))
endef

export CXG_SERVER_PORT := $(call get_or_else_dev_env_default,CXG_SERVER_PORT)
export CXG_CLIENT_PORT := $(call get_or_else_dev_env_default,CXG_CLIENT_PORT)
export JEST_ENV := $(call get_or_else_dev_env_default,JEST_ENV)
export DATASET := $(call get_or_else_dev_env_default,DATASET)
export CXG_OPTIONS := $(call get_or_else_dev_env_default,CXG_OPTIONS)

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
