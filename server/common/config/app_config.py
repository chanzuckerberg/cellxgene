import yaml
from flatten_dict import unflatten

from server import display_version as cellxgene_display_version
from server.common.config.dataset_config import DatasetConfig
from server.common.config.server_config import ServerConfig
from server.common.default_config import get_default_config
from server.common.errors import ConfigurationError



class AppConfig(object):
    """AppConfig stores all the configuration for cellxgene.  The configuration is divided into two main parts:
    server attributes, and dataset attributes. The server_config contains attributes that refer to the server process
    as a whole.   The default_dataset_config referes to attributes that are associated with the features and
    presentations of a dataset.  The dataset config attributes can be overridden depending on the url by which the
    dataset was accessed.  These are stored in dataroot_config.
    AppConfig has methods to initialize, modify, and access the configuration.
    """

    def __init__(self):

        # the default configuration (see default_config.py)
        self.default_config = get_default_config()
        # the server configuration
        self.server_config = ServerConfig(self, self.default_config["server"])
        # the dataset config, unless overridden by an entry in dataroot_config
        self.default_dataset_config = DatasetConfig(None, self, self.default_config["dataset"])
        # a dictionary of keys to DatasetConfig objects.  Each key must exist in the multi_dataset__dataroot
        # attribute of the server_config.
        self.dataroot_config = {}

        # Set to true when config_completed is called
        self.is_completed = False

    def get_dataset_config(self, dataroot_key):
        if self.server_config.single_dataset__datapath:
            return self.default_dataset_config
        else:
            return self.dataroot_config.get(dataroot_key, self.default_dataset_config)

    def check_config(self):
        """Verify all the attributes have been checked"""
        if not self.is_completed:
            raise ConfigurationError("The configuration has not been completed")
        self.server_config.check_config()
        self.default_dataset_config.check_config()
        for dataset_config in self.dataroot_config.values():
            dataset_config.check_config()

    def update_server_config(self, **kw):
        self.server_config.update(**kw)
        self.is_complete = False

    def update_default_dataset_config(self, **kw):
        self.default_dataset_config.update(**kw)
        # update all the other dataset configs, if any
        for value in self.dataroot_config.values():
            value.update(**kw)
        self.is_complete = False

    def update_from_config_file(self, config_file):
        with open(config_file) as fyaml:
            config = yaml.load(fyaml, Loader=yaml.FullLoader)

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

        self.is_complete = False

    def write_config(self, config_file):
        """output the config to a yaml file"""
        server = self.server_config.create_mapping(self.server_config.default_config)
        dataset = self.default_dataset_config.create_mapping(self.default_dataset_config.default_config)
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

        config = unflatten(config, splitter=lambda key: key.split("__"))
        yaml.dump(config, open(config_file, "w"))

    def changes_from_default(self):
        """Return all the attribute that are different from the default"""
        diff_server = self.server_config.changes_from_default()
        diff_dataset = self.default_dataset_config.changes_from_default()
        diff = dict(server=diff_server, dataset=diff_dataset)
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

    def get_client_config(self, data_adaptor):
        """
        Return the configuration as required by the /config REST route
        """

        server_config = self.server_config
        dataset_config = data_adaptor.dataset_config
        annotation = dataset_config.user_annotations
        auth = server_config.auth

        # FIXME The current set of config is not consistently presented:
        # we have camalCase, hyphen-text, and underscore_text

        # make sure the configuration has been checked.
        self.check_config()

        # features
        features = [f.todict() for f in data_adaptor.get_features(annotation)]

        # display_names
        title = self.get_title(data_adaptor)
        about = self.get_about(data_adaptor)

        display_names = dict(engine=data_adaptor.get_name(), dataset=title)

        # library_versions
        library_versions = {}
        library_versions.update(data_adaptor.get_library_versions())
        library_versions["cellxgene"] = cellxgene_display_version

        # links
        links = {"about-dataset": about}

        # parameters
        parameters = {
            "layout": dataset_config.embeddings__names,
            "max-category-items": dataset_config.presentation__max_categories,
            "obs_names": server_config.single_dataset__obs_names,
            "var_names": server_config.single_dataset__var_names,
            "diffexp_lfc_cutoff": dataset_config.diffexp__lfc_cutoff,
            "backed": server_config.adaptor__anndata_adaptor__backed,
            "disable-diffexp": not dataset_config.diffexp__enable,
            "enable-reembedding": dataset_config.embeddings__enable_reembedding,
            "annotations": False,
            "annotations_file": None,
            "annotations_dir": None,
            "annotations_cell_ontology_enabled": False,
            "annotations_cell_ontology_obopath": None,
            "annotations_cell_ontology_terms": None,
            "custom_colors": dataset_config.presentation__custom_colors,
            "diffexp-may-be-slow": False,
            "about_legal_tos": dataset_config.app__about_legal_tos,
            "about_legal_privacy": dataset_config.app__about_legal_privacy,
        }

        # corpora dataset_props
        # TODO/Note: putting info from the dataset into the /config is not ideal.
        # However, it is definitely not part of /schema, and we do not have a top-level
        # route for data properties.  Consider creating one at some point.
        corpora_props = data_adaptor.get_corpora_props()
        if corpora_props and "default_embedding" in corpora_props:
            default_embedding = corpora_props["default_embedding"]
            if isinstance(default_embedding, str) and default_embedding.startswith("X_"):
                default_embedding = default_embedding[2:]  # drop X_ prefix
            if default_embedding in data_adaptor.get_embedding_names():
                parameters["default_embedding"] = default_embedding

        data_adaptor.update_parameters(parameters)
        if annotation:
            annotation.update_parameters(parameters, data_adaptor)

        # gather it all together
        c = {}
        config = c["config"] = {}
        config["features"] = features
        config["displayNames"] = display_names
        config["library_versions"] = library_versions
        config["links"] = links
        config["parameters"] = parameters
        config["corpora_props"] = corpora_props
        config["limits"] = {
            "column_request_max": server_config.limits__column_request_max,
            "diffexp_cellcount_max": server_config.limits__diffexp_cellcount_max,
        }

        if dataset_config.app__authentication_enable and auth.is_valid_authentication_type():
            config["authentication"] = {
                "requires_client_login": auth.requires_client_login(),
            }
            if auth.requires_client_login():
                config["authentication"].update({
                    "login": auth.get_login_url(data_adaptor),
                    "logout": auth.get_logout_url(data_adaptor),
                })

        return c

    def get_client_userinfo(self, data_adaptor):
        """
        Return the userinfo as required by the /userinfo REST route
        """

        server_config = self.server_config
        dataset_config = data_adaptor.dataset_config
        auth = server_config.auth

        # make sure the configuration has been checked.
        self.check_config()

        if dataset_config.app__authentication_enable and auth.is_valid_authentication_type():
            userinfo = {}
            userinfo["userinfo"] = {
                "is_authenticated": auth.is_user_authenticated(),
                "username": auth.get_user_name(),
                "user_id": auth.get_user_id()
            }
            return userinfo
        else:
            return None
