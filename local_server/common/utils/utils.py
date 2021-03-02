import importlib.util
import logging
import pkgutil

from local_server.common.errors import ConfigurationError


def import_plugins(plugin_module):
    """
    Load optional plugin modules from local_server.common.plugins

    If you would like to customize cellxgene, you can add submodules to server.common.plugins before running the app.
    This code will import each, loading the code in each. If no plugins are defined, initializing the app continues as
    normal.
    """
    loaded_modules = []
    try:
        pkg = importlib.import_module(plugin_module)
        for loader, name, is_pkg in pkgutil.walk_packages(pkg.__path__):
            full_name = f"{plugin_module}.{name}"
            try:
                module = importlib.import_module(full_name)
            except Exception as e:
                raise ConfigurationError(f"Unexpected error while importing plugin: {plugin_module}.{name}: {str(e)}")
            loaded_modules.append(module)
    except ModuleNotFoundError as e:
        #  This exception occurs when the plugin_module does not exist (not an error).
        logging.debug(f"No plugins found in module: {plugin_module}: {str(e)}")

    return loaded_modules
