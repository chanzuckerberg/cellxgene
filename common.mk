PROJECT_ROOT := $(shell git rev-parse --show-toplevel)
PATH := $(PATH):$(PROJECT_ROOT)/scripts

SHELL := env PATH='$(PATH)' /bin/bash

ifeq ($(shell which jq),)
$(error Please install jq using "apt-get install jq" or "brew install jq")
endif

define env_or_else_default
$(if $($(1)),$($(1)),$(shell jq -r '.$(1)' "$(PROJECT_ROOT)/environment.default.json"))
endef

# if not a full path, create a full path relative to the project root
define full_path
$(shell [[ $(1) = /* ]] && echo $(1) || echo $(PROJECT_ROOT)/$(1))
endef

export CXG_SERVER_PORT := $(call env_or_else_default,CXG_SERVER_PORT)
export CXG_CLIENT_PORT := $(call env_or_else_default,CXG_CLIENT_PORT)
export CXG_OPTIONS := $(call env_or_else_default,CXG_OPTIONS)
export DATASET := $(call full_path,$(call env_or_else_default,DATASET))
export JEST_ENV := $(call env_or_else_default,JEST_ENV)

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

