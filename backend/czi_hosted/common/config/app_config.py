import yaml
from flatten_dict import unflatten

from backend.czi_hosted.common.config.external_config import ExternalConfig
from backend.czi_hosted.common.config.dataset_config import DatasetConfig
from backend.czi_hosted.common.config.server_config import ServerConfig
from backend.czi_hosted.common.errors import ConfigurationError
from backend.czi_hosted.default_config import get_default_config


class AppConfig(object):
    """
    AppConfig stores all the configuration for cellxgene.
    AppConfig contains one or more DatasetConfig(s) and one ServerConfig.
    The server_config contains attributes that refer to the server process as a whole.
    The default_dataset_config refers to attributes that are associated with the features and
    presentations of a dataset.
    The dataset config attributes can be overridden depending on the url by which the
    dataset was accessed.  These are stored in dataroot_config.
    AppConfig has methods to initialize, modify, and access the configuration.
    """

    def __init__(self):

        # the default configuration (see default_config.py)
        # TODO @madison -- if we always read from the default config (hard coded path) can we set those values as
        #  defaults within the config class?
        self.default_config = get_default_config()
        # the server configuration
        self.server_config = ServerConfig(self, self.default_config["server"])
        # the dataset config, unless overridden by an entry in dataroot_config
        self.default_dataset_config = DatasetConfig(None, self, self.default_config["dataset"])
        # a dictionary of keys to DatasetConfig objects.  Each key must exist in the multi_dataset__dataroot
        # attribute of the server_config. The default dataset config will apply to all datasets unless a different set
        # of config vars was passed for a specific dataset under the multidataset config. For example:
        """
        per_dataset_config:
            d1:
               user_annotations:
                   enable:  false
            d2:
               user_annotations:
                   enable:  true
         """
        #   dataroot config
        self.dataroot_config = {}

        #  external config
        self.external_config = ExternalConfig(self, self.default_config["external"])

        # Set to true when config_completed is called
        self.is_completed = False

    def get_dataset_config(self, dataroot_key):
        if self.server_config.single_dataset__datapath:
            return self.default_dataset_config
        else:
            return self.dataroot_config.get(dataroot_key, self.default_dataset_config)

    def check_config(self):
        """Verify all the attributes in the config have been type checked"""
        if not self.is_completed:
            raise ConfigurationError("The configuration has not been completed")
        self.server_config.check_config()
        self.default_dataset_config.check_config()
        for dataset_config in self.dataroot_config.values():
            dataset_config.check_config()
        self.external_config.check_config()

    def update_server_config(self, **kw):
        self.server_config.update(**kw)
        self.is_complete = False

    def update_default_dataset_config(self, **kw):
        self.default_dataset_config.update(**kw)
        # update all the other dataset configs, if any
        for value in self.dataroot_config.values():
            value.update(**kw)
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

        if len(path) < 1 or path[0] not in ("server", "dataset", "per_dataset_config"):
            raise ConfigurationError("path must start with 'server', 'dataset', or 'per_dataset_config'")

        if path[0] == "server":
            attr = "__".join(path[1:])
            try:
                self.update_server_config(**{attr: value})
            except ConfigurationError:
                raise ConfigurationError(f"unknown config parameter at path: '{str(path)}'")
        elif path[0] == "dataset":
            attr = "__".join(path[1:])
            try:
                self.update_default_dataset_config(**{attr: value})
            except ConfigurationError:
                raise ConfigurationError(f"unknown config parameter at path: '{str(path)}'")

        elif path[0] == "per_dataset_config":
            if len(path) < 2:
                raise ConfigurationError(f"missing dataroot when using per_dataset_config: got '{path}'")
            dataroot = path[1]
            if dataroot not in self.dataroot_config:
                dataroots = str(list(self.dataroot_config.keys()))
                raise ConfigurationError(
                    f"unknown dataroot when using per_dataset_config: got '{path}',"
                    f" dataroots specified in config are {dataroots}"
                )

            attr = "__".join(path[2:])
            try:
                self.dataroot_config[dataroot].update(**{attr: value})
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
            self.default_dataset_config.update_from_config(config["dataset"], "dataset")

        per_dataset_config = config.get("per_dataset_config", {})
        for key, dataroot_config in per_dataset_config.items():
            # first create and initialize the dataroot with the default config
            self.add_dataroot_config(key, **config["dataset"])
            # then apply the per dataset configuration
            self.dataroot_config[key].update_from_config(dataroot_config, f"per_dataset_config__{key}")

        if config.get("external"):
            self.external_config.update_from_config(config["external"], "external")

        self.is_complete = False

    def config_to_dict(self):
        """return the configuration as an unflattened dict"""
        server = self.server_config.create_mapping(self.server_config.default_config)
        dataset = self.default_dataset_config.create_mapping(self.default_dataset_config.default_config)
        external = self.external_config.create_mapping(self.external_config.default_config)
        config = dict(server={}, dataset={})
        for attrname in server.keys():
            config["server__" + attrname] = getattr(self.server_config, attrname)
        for attrname in dataset.keys():
            config["dataset__" + attrname] = getattr(self.default_dataset_config, attrname)
        if self.dataroot_config:
            config["per_dataset_config"] = {}
        for dataroot_tag, dataroot_config in self.dataroot_config.items():
            dataset = dataroot_config.create_mapping(dataroot_config.default_config)
            for attrname in dataset.keys():
                config[f"per_dataset_config__{dataroot_tag}__" + attrname] = getattr(dataroot_config, attrname)
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
        diff_dataset = self.default_dataset_config.changes_from_default()
        diff_external = self.external.changes_from_default()
        diff = dict(server=diff_server, dataset=diff_dataset, external=diff_external)
        return diff

    def add_dataroot_config(self, dataroot_tag, **kw):
        """Create a new dataset config object based on the default dataset config, and kw parameters"""
        if dataroot_tag in self.dataroot_config:
            raise ConfigurationError(f"dataroot config already exists: {dataroot_tag}")
        if type(self.server_config.multi_dataset__dataroot) != dict:
            raise ConfigurationError("The server__multi_dataset__dataroot must be a dictionary")
        if dataroot_tag not in self.server_config.multi_dataset__dataroot:
            raise ConfigurationError(f"The dataroot_tag ({dataroot_tag}) not found in server__multi_dataset__dataroot")

        self.is_completed = False
        self.dataroot_config[dataroot_tag] = DatasetConfig(dataroot_tag, self, self.default_config["dataset"])
        flat_config = self.default_dataset_config.create_mapping(self.default_dataset_config.default_config)
        config = {key: value[1] for key, value in flat_config.items()}
        self.dataroot_config[dataroot_tag].update(**config)
        self.dataroot_config[dataroot_tag].update_from_config(kw, dataroot_tag)

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
        self.default_dataset_config.complete_config(context)
        for dataroot_config in self.dataroot_config.values():
            dataroot_config.complete_config(context)

        self.is_completed = True
        self.check_config()

    def get_matrix_data_cache_manager(self):
        return self.server_config.matrix_data_cache_manager

    def is_multi_dataset(self):
        return self.server_config.multi_dataset__dataroot is not None

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
