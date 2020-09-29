import logging
import sys
import os
from server.common.utils.utils import import_plugins

__version__ = "0.16.0"
display_version = "cellxgene v" + __version__

try:
    import_plugins("server.plugins")
except Exception as e:
    # Make sure to exit in this case, as the server may not be configured as expected.
    logging.critical(f"Error in import_plugins: {str(e)}")
    sys.exit(1)

try:
    PROJECT_ROOT = os.popen("git rev-parse --show-toplevel").read().strip()
    os.environ["PROJECT_ROOT"] = PROJECT_ROOT
except Exception as e:
    logging.critical(f"Issue setting project root: {str(e)}")
    sys.exit(1)
