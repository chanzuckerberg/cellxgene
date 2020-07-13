from server.common.utils import import_plugins

__version__ = "0.15.0"
display_version = "cellxgene v" + __version__

import_plugins("server.plugins")
