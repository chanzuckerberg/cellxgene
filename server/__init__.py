from server.common.utils import import_plugins

__version__ = "0.15.0"


import_plugins("server.plugins")
