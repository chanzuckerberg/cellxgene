import yaml
from flatten_dict import unflatten

from backend.server.default_config import get_default_config
from backend.server.common.config.dataset_config import DatasetConfig
from backend.server.common.config.server_config import ServerConfig
from backend.server.common.config.external_config import ExternalConfig
from backend.common.errors import ConfigurationError


class AppConfig(object):
    """
    AppConfig stores all the configuration for cellxgene.
    AppConfig contains one or more DatasetConfig(s) and one ServerConfig.
    The server_config contains attributes that refer to the server process as a whole.
    The dataset_config refers to attributes that are associated with the features and
    presentations of a dataset.
    AppConfig has methods to initialize, modify, and access the configuration.
    """

    def __init__(self):

        # the default configuration (see default_config.py)
        # TODO @madison -- if we always read from the default config (hard coded path) can we set those values as
        #  defaults within the config class?
        self.default_config = get_default_config()
        # the server configuration
        self.server_config = ServerConfig(self, self.default_config["server"])
        # the dataset config
        self.dataset_config = DatasetConfig(None, self, self.default_config["dataset"])
        #  external config
        self.external_config = ExternalConfig(self, self.default_config["external"])

        # Set to true when config_completed is called
        self.is_completed = False

    def get_dataset_config(self):
        return self.dataset_config

    def check_config(self):
        """Verify all the attributes in the config have been type checked"""
        if not self.is_completed:
            raise ConfigurationError("The configuration has not been completed")
        self.server_config.check_config()
        self.dataset_config.check_config()
        self.external_config.check_config()

    def update_server_config(self, **kw):
        self.server_config.update(**kw)
        self.is_complete = False

    def update_dataset_config(self, **kw):
        self.dataset_config.update(**kw)
        self.is_complete = False

    def update_single_config_from_path_and_value(self, path, value):
        """Update a single config parameter with the value.
        Path is a list of string, that gives a path to the config parameter to be updated.
        For example, path may be ["server","app","port"].
        """
        self.is_complete = False
        if not isinstance(path, list):
            raise ConfigurationError(f"path must be a list of strings, got '{str(path)}'")
        for part in path:
            if not isinstance(part, str):
                raise ConfigurationError(f"path must be a list of strings, got '{str(path)}'")

        if len(path) < 1 or path[0] not in ("server", "dataset"):
            raise ConfigurationError("path must start with 'server', or 'dataset'")

        if path[0] == "server":
            attr = "__".join(path[1:])
            try:
                self.update_server_config(**{attr: value})
            except ConfigurationError:
                raise ConfigurationError(f"unknown config parameter at path: '{str(path)}'")
        elif path[0] == "dataset":
            attr = "__".join(path[1:])
            try:
                self.update_dataset_config(**{attr: value})
            except ConfigurationError:
                raise ConfigurationError(f"unknown config parameter at path: '{str(path)}'")

    def update_from_config_file(self, config_file):
        try:
            with open(config_file) as yml_file:
                config = yaml.safe_load(yml_file)
        except yaml.YAMLError as e:
            raise ConfigurationError(f"The specified config file contained an error: {e}")
        except OSError as e:
            raise ConfigurationError(f"Issue retrieving the specified config file: {e}")

        if config.get("server"):
            self.server_config.update_from_config(config["server"], "server")
        if config.get("dataset"):
            self.dataset_config.update_from_config(config["dataset"], "dataset")

        if config.get("external"):
            self.external_config.update_from_config(config["external"], "external")

        self.is_complete = False

    def config_to_dict(self):
        """return the configuration as an unflattened dict"""
        server = self.server_config.create_mapping(self.server_config.default_config)
        dataset = self.dataset_config.create_mapping(self.dataset_config.default_config)
        external = self.external_config.create_mapping(self.external_config.default_config)
        config = dict(server={}, dataset={})
        for attrname in server.keys():
            config["server__" + attrname] = getattr(self.server_config, attrname)
        for attrname in dataset.keys():
            config["dataset__" + attrname] = getattr(self.dataset_config, attrname)
        for attrname in external.keys():
            config["external__" + attrname] = getattr(self.external_config, attrname)

        config = unflatten(config, splitter=lambda key: key.split("__"))
        return config

    def write_config(self, config_file):
        """output the config to a yaml file"""
        config = self.config_to_dict()
        yaml.dump(config, open(config_file, "w"))

    def changes_from_default(self):
        """Return all the attribute that are different from the default"""
        diff_server = self.server_config.changes_from_default()
        diff_dataset = self.dataset_config.changes_from_default()
        diff_external = self.external.changes_from_default()
        diff = dict(server=diff_server, dataset=diff_dataset, external=diff_external)
        return diff

    def complete_config(self, messagefn=None):
        """The configure options are checked, and any additional setup based on the config
        parameters is done"""

        if messagefn is None:

            def noop(message):
                pass

            messagefn = noop

        # TODO: to give better error messages we can add a mapping between where each config
        # attribute originated (e.g. command line argument or config file), then in the error
        # messages we can give correct context for attributes with bad value.
        context = dict(messagefn=messagefn)

        # complete config for external_config first, since this may update values in the other sections
        self.external_config.complete_config(context)
        self.server_config.complete_config(context)
        self.dataset_config.complete_config(context)

        self.is_completed = True
        self.check_config()

    def get_matrix_data_cache_manager(self):
        return self.server_config.matrix_data_cache_manager

    def get_title(self, data_adaptor):
        return (
            self.server_config.single_dataset__title
            if self.server_config.single_dataset__title
            else data_adaptor.get_title()
        )

    def get_about(self, data_adaptor):
        return (
            self.server_config.single_dataset__about
            if self.server_config.single_dataset__about
            else data_adaptor.get_about()
        )
