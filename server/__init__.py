from server.common.utils import import_plugins
import logging
import sys

__version__ = "0.16.0"
display_version = "cellxgene v" + __version__

try:
    import_plugins("server.plugins")
except Exception as e:
    # Make sure to exit in this case, as the server may not be configured as expected.
    logging.critical(f"Error in import_plugins: {str(e)}")
    sys.exit(1)
