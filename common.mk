SHELL=/bin/bash

PROJECT_ROOT := $(shell git rev-parse --show-toplevel)

# get_or_else_dev_env_default
# - If a variable is defined, return its value
# - Else return the default value from environment.dev
define get_or_else_dev_env_default
$(if $($(1)),$($(1)),$(shell sed -n 's/$(1)=\(.*\)/\1/p' < $(PROJECT_ROOT)/environment.dev))
endef

CXG_SERVER_PORT := $(call get_or_else_dev_env_default,CXG_SERVER_PORT)
CXG_CLIENT_PORT := $(call get_or_else_dev_env_default,CXG_CLIENT_PORT)
JEST_ENV := $(call get_or_else_dev_env_default,JEST_ENV)
