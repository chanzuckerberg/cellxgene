import logging
import sys
from server.common.utils.utils import import_plugins

__version__ = "1.0.1"
display_version = "CellxGene v" + __version__

try:
    import_plugins("server.plugins")
except Exception as e:
    # Make sure to exit in this case, as the server may not be configured as expected.
    logging.critical(f"Error in import_plugins: {str(e)}")
    sys.exit(1)
