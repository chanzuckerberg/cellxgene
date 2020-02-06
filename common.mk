SHELL=/bin/bash

# CHECK ENVIRONMENT
ifndef CXG_SERVER_PORT
$(error Please run "source environment.<env>" in the root directory before running make commands)
endif
