import copy
import os
import sys
import warnings
from os.path import splitext, basename, isdir
from urllib.parse import urlparse, quote_plus

import yaml
from flatten_dict import flatten, unflatten

import server.compute.diffexp_cxg as diffexp_tiledb
from server import display_version as cellxgene_display_version
from server.auth.auth import AuthTypeFactory
from server.common.annotations.hosted_tiledb import AnnotationsHostedTileDB
from server.common.annotations.local_file_csv import AnnotationsLocalFile
from server.common.data_locator import discover_s3_region_name
from server.common.default_config import get_default_config
from server.common.errors import ConfigurationError, DatasetAccessError, OntologyLoadFailure
from server.common.utils.utils import custom_format_warning, find_available_port, is_port_available
from server.data_common.matrix_loader import MatrixDataLoader, MatrixDataCacheManager, MatrixDataType
from server.db.db_utils import DbUtils

DEFAULT_SERVER_PORT = 5005
# anything bigger than this will generate a special message
BIG_FILE_SIZE_THRESHOLD = 100 * 2 ** 20  # 100MB


class AppFeature(object):
    def __init__(self, path, available=False, method="POST", extra={}):
        self.path = path
        self.available = available
        self.method = method
        self.extra = extra
        for k, v in extra.items():
            setattr(self, k, v)

    def todict(self):
        d = dict(available=self.available, method=self.method, path=self.path)
        d.update(self.extra)
        return d


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

        self.server_config.update_from_config(config["server"], "server")
        self.default_dataset_config.update_from_config(config["dataset"], "dataset")

        per_dataset_config = config.get("per_dataset_config", {})
        for key, dataroot_config in per_dataset_config.items():
            self.add_dataroot_config(key, **dataroot_config)

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


class BaseConfig(object):
    """This class handles the mechanics of updating and checking attributes.
    Derived classes are expected to store the actual attributes"""

    def __init__(self, app_config, default_config, dictval_cases={}):
        # reference back to the app_config
        self.app_config = app_config
        # the complete set of attribute and their default values (unflattened)
        self.default_config = default_config
        # attributes where the value may be a dict (and therefore are not flattened)
        self.dictval_cases = dictval_cases
        # used to make sure every attribute value is checked
        self.attr_checked = {k: False for k in self.create_mapping(default_config).keys()}

    def create_mapping(self, config):
        """Create a mapping from attribute names to (location in the config tree, value)"""
        dc = copy.deepcopy(config)
        mapping = {}

        # special cases where the value could be a dict.
        # If its value is not None, the entry is added to the mapping, and not included
        # in the flattening below.
        for dictval_case in self.dictval_cases:
            cur = dc
            for part in dictval_case[:-1]:
                cur = cur.get(part, {})
            val = cur.get(dictval_case[-1])
            if val is not None:
                key = "__".join(dictval_case)
                mapping[key] = (dictval_case, val)
                del cur[dictval_case[-1]]

        flat_config = flatten(dc)
        for key, value in flat_config.items():
            # name of the attribute
            attr = "__".join(key)
            mapping[attr] = (key, value)

        return mapping

    def check_attr(self, attrname, vtype):
        val = getattr(self, attrname)
        if type(vtype) in (list, tuple):
            if type(val) not in vtype:
                tnames = ",".join([x.__name__ for x in vtype])
                raise ConfigurationError(
                    f"Invalid type for attribute: {attrname}, expected types ({tnames}), got {type(val).__name__}"
                )
        else:
            if type(val) != vtype:
                raise ConfigurationError(
                    f"Invalid type for attribute: {attrname}, "
                    f"expected type {vtype.__name__}, got {type(val).__name__}"
                )

        self.attr_checked[attrname] = True

    def check_config(self):
        mapping = self.create_mapping(self.default_config)
        for key in mapping.keys():
            if not self.attr_checked[key]:
                raise ConfigurationError(f"The attr '{key}' has not been checked")

    def update(self, **kw):
        for key, value in kw.items():
            if not hasattr(self, key):
                raise ConfigurationError(f"unknown config parameter {key}.")
            try:
                if type(value) == tuple:
                    # convert tuple values to list values
                    value = list(value)
                setattr(self, key, value)
            except KeyError:
                raise ConfigurationError(f"Unable to set config parameter {key}.")

            self.attr_checked[key] = False

    def update_from_config(self, config, prefix):
        mapping = self.create_mapping(config)
        for attr, (key, value) in mapping.items():
            if not hasattr(self, attr):
                raise ConfigurationError(f"Unknown key from config file: {prefix}__{attr}")
            try:
                setattr(self, attr, value)
            except KeyError:
                raise ConfigurationError(f"Unable to set config attribute: {prefix}__{attr}")

            self.attr_checked[attr] = False

    def changes_from_default(self):
        """Return all the attribute that are different from the default"""
        mapping = self.create_mapping(self.default_config)
        diff = []
        for attrname, (key, defval) in mapping.items():
            curval = getattr(self, attrname)
            if curval != defval:
                diff.append((attrname, curval, defval))
        return diff


class ServerConfig(BaseConfig):
    """Manages the config attribute associated with the server."""

    def __init__(self, app_config, default_config):
        dictval_cases = [
            ("app", "csp_directives"),
            ("authentication", "params_oauth", "cookie"),
            ("adaptor", "cxg_adaptor", "tiledb_ctx"),
            ("multi_dataset", "dataroot"),
        ]
        super().__init__(app_config, default_config, dictval_cases)

        dc = default_config
        try:
            self.app__verbose = dc["app"]["verbose"]
            self.app__debug = dc["app"]["debug"]
            self.app__host = dc["app"]["host"]
            self.app__port = dc["app"]["port"]
            self.app__open_browser = dc["app"]["open_browser"]
            self.app__force_https = dc["app"]["force_https"]
            self.app__flask_secret_key = dc["app"]["flask_secret_key"]
            self.app__generate_cache_control_headers = dc["app"]["generate_cache_control_headers"]
            self.app__server_timing_headers = dc["app"]["server_timing_headers"]
            self.app__csp_directives = dc["app"]["csp_directives"]

            self.authentication__type = dc["authentication"]["type"]
            self.authentication__params_oauth__api_base_url = dc["authentication"]["params_oauth"]["api_base_url"]
            self.authentication__params_oauth__client_id = dc["authentication"]["params_oauth"]["client_id"]
            self.authentication__params_oauth__client_secret = dc["authentication"]["params_oauth"]["client_secret"]
            self.authentication__params_oauth__callback_base_url = \
                dc["authentication"]["params_oauth"]["callback_base_url"]
            self.authentication__params_oauth__session_cookie = dc["authentication"]["params_oauth"]["session_cookie"]
            self.authentication__params_oauth__cookie = dc["authentication"]["params_oauth"]["cookie"]

            self.multi_dataset__dataroot = dc["multi_dataset"]["dataroot"]
            self.multi_dataset__index = dc["multi_dataset"]["index"]
            self.multi_dataset__allowed_matrix_types = dc["multi_dataset"]["allowed_matrix_types"]
            self.multi_dataset__matrix_cache__max_datasets = dc["multi_dataset"]["matrix_cache"]["max_datasets"]
            self.multi_dataset__matrix_cache__timelimit_s = dc["multi_dataset"]["matrix_cache"]["timelimit_s"]

            self.single_dataset__datapath = dc["single_dataset"]["datapath"]
            self.single_dataset__obs_names = dc["single_dataset"]["obs_names"]
            self.single_dataset__var_names = dc["single_dataset"]["var_names"]
            self.single_dataset__about = dc["single_dataset"]["about"]
            self.single_dataset__title = dc["single_dataset"]["title"]

            self.diffexp__alg_cxg__max_workers = dc["diffexp"]["alg_cxg"]["max_workers"]
            self.diffexp__alg_cxg__cpu_multiplier = dc["diffexp"]["alg_cxg"]["cpu_multiplier"]
            self.diffexp__alg_cxg__target_workunit = dc["diffexp"]["alg_cxg"]["target_workunit"]

            self.data_locator__s3__region_name = dc["data_locator"]["s3"]["region_name"]

            self.adaptor__cxg_adaptor__tiledb_ctx = dc["adaptor"]["cxg_adaptor"]["tiledb_ctx"]
            self.adaptor__anndata_adaptor__backed = dc["adaptor"]["anndata_adaptor"]["backed"]

            self.limits__diffexp_cellcount_max = dc["limits"]["diffexp_cellcount_max"]
            self.limits__column_request_max = dc["limits"]["column_request_max"]

        except KeyError as e:
            raise ConfigurationError(f"Unexpected config: {str(e)}")

        # The matrix data cache manager is created during the complete_config and stored here.
        self.matrix_data_cache_manager = None

        # The authentication object
        self.auth = None

    def complete_config(self, context):
        self.handle_app(context)
        self.handle_data_source(context)
        self.handle_authentication(context)
        self.handle_data_locator(context)
        self.handle_adaptor(context)  # may depend on data_locator
        self.handle_single_dataset(context)  # may depend on adaptor
        self.handle_multi_dataset(context)  # may depend on adaptor
        self.handle_diffexp(context)
        self.handle_limits(context)

        self.check_config()

    def handle_app(self, context):
        self.check_attr("app__verbose", bool)
        self.check_attr("app__debug", bool)
        self.check_attr("app__host", str)
        self.check_attr("app__port", (type(None), int))
        self.check_attr("app__open_browser", bool)
        self.check_attr("app__force_https", bool)
        self.check_attr("app__flask_secret_key", (type(None), str))
        self.check_attr("app__generate_cache_control_headers", bool)
        self.check_attr("app__server_timing_headers", bool)
        self.check_attr("app__csp_directives", (type(None), dict))

        if self.app__port:
            try:
                if not is_port_available(self.app__host, self.app__port):
                    raise ConfigurationError(
                        f"The port selected {self.app__port} is in use, please configure an open port."
                    )
            except OverflowError:
                raise ConfigurationError(f"Invalid port: {self.app__port}")
        else:
            try:
                default_server_port = int(os.environ.get("CXG_SERVER_PORT", DEFAULT_SERVER_PORT))
            except ValueError:
                raise ConfigurationError(
                    "Invalid port from environment variable CXG_SERVER_PORT: " + os.environ.get("CXG_SERVER_PORT")
                )
            try:
                self.app__port = find_available_port(self.app__host, default_server_port)
            except OverflowError:
                raise ConfigurationError(f"Invalid port: {default_server_port}")

        if self.app__debug:
            context["messagefn"]("in debug mode, setting verbose=True and open_browser=False")
            self.app__verbose = True
            self.app__open_browser = False
        else:
            warnings.formatwarning = custom_format_warning

        if not self.app__verbose:
            sys.tracebacklimit = 0

        # secret key:
        #   first, from CXG_SECRET_KEY environment variable
        #   second, from config file
        self.app__flask_secret_key = os.environ.get("CXG_SECRET_KEY", self.app__flask_secret_key)

        # CSP Directives are a dict of string: list(string) or string: string
        if self.app__csp_directives is not None:
            for k, v in self.app__csp_directives.items():
                if not isinstance(k, str):
                    raise ConfigurationError("CSP directive names must be a string.")
                if isinstance(v, list):
                    for policy in v:
                        if not isinstance(policy, str):
                            raise ConfigurationError("CSP directive value must be a string or list of strings.")
                elif not isinstance(v, str):
                    raise ConfigurationError("CSP directive value must be a string or list of strings.")

    def handle_authentication(self, context):
        self.check_attr("authentication__type", (type(None), str))

        # oauth
        ptypes = str if self.authentication__type == "oauth" else (type(None), str)
        self.check_attr("authentication__params_oauth__api_base_url", ptypes)
        self.check_attr("authentication__params_oauth__client_id", ptypes)
        self.check_attr("authentication__params_oauth__client_secret", ptypes)
        self.check_attr("authentication__params_oauth__callback_base_url", (type(None), str))
        self.check_attr("authentication__params_oauth__session_cookie", bool)

        if self.authentication__params_oauth__session_cookie:
            self.check_attr("authentication__params_oauth__cookie", (type(None), dict))
        else:
            self.check_attr("authentication__params_oauth__cookie", dict)
        #   secret key: first, from CXG_OAUTH_CLIENT_SECRET environment variable
        #   second, from config file
        self.authentication__params__oauth__client_secret = os.environ.get(
            "CXG_OAUTH_CLIENT_SECRET", self.authentication__params_oauth__client_secret)

        self.auth = AuthTypeFactory.create(self.authentication__type, self)
        if self.auth is None:
            raise ConfigurationError(f"Unknown authentication type: {self.authentication__type}")

    def handle_data_locator(self, context):
        self.check_attr("data_locator__s3__region_name", (type(None), bool, str))
        if self.data_locator__s3__region_name is True:
            path = self.single_dataset__datapath or self.multi_dataset__dataroot
            if type(path) == dict:
                # if multi_dataset__dataroot is a dict, then use the first key
                # that is in s3.   NOTE:  it is not supported to have dataroots
                # in different regions.
                paths = [val.get("dataroot") for val in path.values()]
                for path in paths:
                    if path.startswith("s3://"):
                        break
            if path.startswith("s3://"):
                region_name = discover_s3_region_name(path)
                if region_name is None:
                    raise ConfigurationError(f"Unable to discover s3 region name from {path}")
            else:
                region_name = None
            self.data_locator__s3__region_name = region_name

    def handle_data_source(self, context):
        self.check_attr("single_dataset__datapath", (str, type(None)))
        self.check_attr("multi_dataset__dataroot", (type(None), dict, str))

        if self.single_dataset__datapath is None:
            if self.multi_dataset__dataroot is None:
                # TODO:  change the error message once dataroot is fully supported
                raise ConfigurationError("missing datapath")
            return
        else:
            if self.multi_dataset__dataroot is not None:
                raise ConfigurationError("must supply only one of datapath or dataroot")

    def handle_single_dataset(self, context):
        self.check_attr("single_dataset__datapath", (str, type(None)))
        self.check_attr("single_dataset__title", (str, type(None)))
        self.check_attr("single_dataset__about", (str, type(None)))
        self.check_attr("single_dataset__obs_names", (str, type(None)))
        self.check_attr("single_dataset__var_names", (str, type(None)))

        if self.single_dataset__datapath is None:
            return

        # create the matrix data cache manager:
        if self.matrix_data_cache_manager is None:
            self.matrix_data_cache_manager = MatrixDataCacheManager(max_cached=1, timelimit_s=None)

        # preload this data set
        matrix_data_loader = MatrixDataLoader(self.single_dataset__datapath, app_config=self.app_config)
        try:
            matrix_data_loader.pre_load_validation()
        except DatasetAccessError as e:
            raise ConfigurationError(str(e))

        file_size = matrix_data_loader.file_size()
        file_basename = basename(self.single_dataset__datapath)
        if file_size > BIG_FILE_SIZE_THRESHOLD:
            context["messagefn"](f"Loading data from {file_basename}, this may take a while...")
        else:
            context["messagefn"](f"Loading data from {file_basename}.")

        if self.single_dataset__about:

            def url_check(url):
                try:
                    result = urlparse(url)
                    if all([result.scheme, result.netloc]):
                        return True
                    else:
                        return False
                except ValueError:
                    return False

            if not url_check(self.single_dataset__about):
                raise ConfigurationError(
                    "Must provide an absolute URL for --about. (Example format: http://example.com)"
                )

    def handle_multi_dataset(self, context):
        self.check_attr("multi_dataset__dataroot", (type(None), dict, str))
        self.check_attr("multi_dataset__index", (type(None), bool, str))
        self.check_attr("multi_dataset__allowed_matrix_types", list)
        self.check_attr("multi_dataset__matrix_cache__max_datasets", int)
        self.check_attr("multi_dataset__matrix_cache__timelimit_s", (type(None), int, float))

        if self.multi_dataset__dataroot is None:
            return

        if type(self.multi_dataset__dataroot) == str:
            default_dict = dict(base_url="d", dataroot=self.multi_dataset__dataroot)
            self.multi_dataset__dataroot = dict(d=default_dict)

        for tag, dataroot_dict in self.multi_dataset__dataroot.items():
            if "base_url" not in dataroot_dict:
                raise ConfigurationError(f"error in multi_dataset__dataroot: missing base_url for tag {tag}")
            if "dataroot" not in dataroot_dict:
                raise ConfigurationError(f"error in multi_dataset__dataroot: missing dataroot, for tag {tag}")

            base_url = dataroot_dict["base_url"]

            # sanity check for well formed base urls
            bad = False
            if type(base_url) != str:
                bad = True
            elif os.path.normpath(base_url) != base_url:
                bad = True
            else:
                base_url_parts = base_url.split("/")
                if [quote_plus(part) for part in base_url_parts] != base_url_parts:
                    bad = True
                if ".." in base_url_parts:
                    bad = True
            if bad:
                raise ConfigurationError(f"error in multi_dataset__dataroot base_url {base_url} for tag {tag}")

        # verify all the base_urls are unique
        base_urls = [d["base_url"] for d in self.multi_dataset__dataroot.values()]
        if len(base_urls) > len(set(base_urls)):
            raise ConfigurationError("error in multi_dataset__dataroot:  base_urls must be unique")

        # error checking
        for mtype in self.multi_dataset__allowed_matrix_types:
            try:
                MatrixDataType(mtype)
            except ValueError:
                raise ConfigurationError(f'Invalid matrix type in "allowed_matrix_types": {mtype}')

        # create the matrix data cache manager:
        if self.matrix_data_cache_manager is None:
            self.matrix_data_cache_manager = MatrixDataCacheManager(
                max_cached=self.multi_dataset__matrix_cache__max_datasets,
                timelimit_s=self.multi_dataset__matrix_cache__timelimit_s,
            )

    def handle_diffexp(self, context):
        self.check_attr("diffexp__alg_cxg__max_workers", (str, int))
        self.check_attr("diffexp__alg_cxg__cpu_multiplier", int)
        self.check_attr("diffexp__alg_cxg__target_workunit", int)

        max_workers = self.diffexp__alg_cxg__max_workers
        cpu_multiplier = self.diffexp__alg_cxg__cpu_multiplier
        cpu_count = os.cpu_count()
        max_workers = min(max_workers, cpu_multiplier * cpu_count)
        diffexp_tiledb.set_config(max_workers, self.diffexp__alg_cxg__target_workunit)

    def handle_adaptor(self, context):
        # cxg
        self.check_attr("adaptor__cxg_adaptor__tiledb_ctx", dict)
        regionkey = "vfs.s3.region"
        if regionkey not in self.adaptor__cxg_adaptor__tiledb_ctx:
            if type(self.data_locator__s3__region_name) == str:
                self.adaptor__cxg_adaptor__tiledb_ctx[regionkey] = self.data_locator__s3__region_name

        from server.data_cxg.cxg_adaptor import CxgAdaptor

        CxgAdaptor.set_tiledb_context(self.adaptor__cxg_adaptor__tiledb_ctx)

        # anndata
        self.check_attr("adaptor__anndata_adaptor__backed", bool)

    def handle_limits(self, context):
        self.check_attr("limits__diffexp_cellcount_max", (type(None), int))
        self.check_attr("limits__column_request_max", (type(None), int))

    def exceeds_limit(self, limit_name, value):
        limit_value = getattr(self, "limits__" + limit_name, None)
        if limit_value is None:  # disabled
            return False
        return value > limit_value


class DatasetConfig(BaseConfig):
    """Manages the config attribute associated with a dataset."""

    def __init__(self, tag, app_config, default_config):
        super().__init__(app_config, default_config)
        self.tag = tag
        dc = default_config
        try:
            self.app__scripts = dc["app"]["scripts"]
            self.app__inline_scripts = dc["app"]["inline_scripts"]
            self.app__about_legal_tos = dc["app"]["about_legal_tos"]
            self.app__about_legal_privacy = dc["app"]["about_legal_privacy"]
            self.app__authentication_enable = dc["app"]["authentication_enable"]

            self.presentation__max_categories = dc["presentation"]["max_categories"]
            self.presentation__custom_colors = dc["presentation"]["custom_colors"]

            self.user_annotations__enable = dc["user_annotations"]["enable"]
            self.user_annotations__type = dc["user_annotations"]["type"]
            self.user_annotations__local_file_csv__directory = dc["user_annotations"]["local_file_csv"]["directory"]
            self.user_annotations__local_file_csv__file = dc["user_annotations"]["local_file_csv"]["file"]
            self.user_annotations__ontology__enable = dc["user_annotations"]["ontology"]["enable"]
            self.user_annotations__ontology__obo_location = dc["user_annotations"]["ontology"]["obo_location"]
            self.user_annotations__hosted_tiledb_array__db_uri = dc["user_annotations"]["hosted_tiledb_array"]["db_uri"]
            self.user_annotations__hosted_tiledb_array__hosted_file_directory = \
                dc["user_annotations"]["hosted_tiledb_array"]["hosted_file_directory"]  # noqa E501

            self.embeddings__names = dc["embeddings"]["names"]
            self.embeddings__enable_reembedding = dc["embeddings"]["enable_reembedding"]

            self.diffexp__enable = dc["diffexp"]["enable"]
            self.diffexp__lfc_cutoff = dc["diffexp"]["lfc_cutoff"]
            self.diffexp__top_n = dc["diffexp"]["top_n"]

        except KeyError as e:
            raise ConfigurationError(f"Unexpected config: {str(e)}")

        # The annotation object is created during complete_config and stored here.
        self.user_annotations = None

    def complete_config(self, context):
        self.handle_app(context)
        self.handle_presentation(context)
        self.handle_user_annotations(context)
        self.handle_embeddings(context)
        self.handle_diffexp(context)

    def handle_app(self, context):
        self.check_attr("app__scripts", list)
        self.check_attr("app__inline_scripts", list)
        self.check_attr("app__about_legal_tos", (type(None), str))
        self.check_attr("app__about_legal_privacy", (type(None), str))
        self.check_attr("app__authentication_enable", bool)

        # scripts can be string (filename) or dict (attributes).   Convert string to dict.
        scripts = []
        for s in self.app__scripts:
            if isinstance(s, str):
                scripts.append({"src": s})
            elif isinstance(s, dict) and isinstance(s["src"], str):
                scripts.append(s)
            else:
                raise ConfigurationError("Scripts must be string or dict")
        self.app__scripts = scripts

    def handle_presentation(self, context):
        self.check_attr("presentation__max_categories", int)
        self.check_attr("presentation__custom_colors", bool)

    def handle_user_annotations(self, context):
        self.check_attr("user_annotations__enable", bool)
        self.check_attr("user_annotations__type", str)
        self.check_attr("user_annotations__local_file_csv__directory", (type(None), str))
        self.check_attr("user_annotations__local_file_csv__file", (type(None), str))
        self.check_attr("user_annotations__ontology__enable", bool)
        self.check_attr("user_annotations__ontology__obo_location", (type(None), str))
        self.check_attr("user_annotations__hosted_tiledb_array__db_uri", (type(None), str))
        self.check_attr("user_annotations__hosted_tiledb_array__hosted_file_directory", (type(None), str))

        if self.user_annotations__enable:
            server_config = self.app_config.server_config
            if not self.app__authentication_enable:
                raise ConfigurationError("user annotations requires authentication to be enabled")
            if not server_config.auth.is_valid_authentication_type():
                auth_type = server_config.authentication__type
                raise ConfigurationError(f"authentication method {auth_type} is not compatible with user annotations")

            # TODO, replace this with a factory pattern once we have more than one way
            # to do annotations.  currently only local_file_csv
            if self.user_annotations__type == "local_file_csv":
                dirname = self.user_annotations__local_file_csv__directory
                filename = self.user_annotations__local_file_csv__file

                if filename is not None and dirname is not None:
                    raise ConfigurationError("'annotations-file' and 'annotations-dir' may not be used together.")

                if filename is not None:
                    lf_name, lf_ext = splitext(filename)
                    if lf_ext and lf_ext != ".csv":
                        raise ConfigurationError(f"annotation file type must be .csv: {filename}")

                if dirname is not None and not isdir(dirname):
                    try:
                        os.mkdir(dirname)
                    except OSError:
                        raise ConfigurationError("Unable to create directory specified by --annotations-dir")

                self.user_annotations = AnnotationsLocalFile(dirname, filename)

                # if the user has specified a fixed label file, go ahead and validate it
                # so that we can remove errors early in the process.
                server_config = self.app_config.server_config
                if server_config.single_dataset__datapath and self.user_annotations__local_file_csv__file:
                    with server_config.matrix_data_cache_manager.data_adaptor(
                            self.tag, server_config.single_dataset__datapath, self.app_config
                    ) as data_adaptor:
                        data_adaptor.check_new_labels(self.user_annotations.read_labels(data_adaptor))

                if self.user_annotations__ontology__enable or self.user_annotations__ontology__obo_location:
                    try:
                        self.user_annotations.load_ontology(self.user_annotations__ontology__obo_location)
                    except OntologyLoadFailure as e:
                        raise ConfigurationError("Unable to load ontology terms\n" + str(e))
            elif self.user_annotations__type == "hosted_tiledb_array":
                self.check_attr("user_annotations__hosted_tiledb_array__db_uri", str)
                self.check_attr("user_annotations__hosted_tiledb_array__hosted_file_directory", str)
                self.user_annotations = AnnotationsHostedTileDB(
                    directory_path=self.user_annotations__hosted_tiledb_array__hosted_file_directory,
                    db=DbUtils(self.user_annotations__hosted_tiledb_array__db_uri),
                )
            else:
                raise ConfigurationError('The only annotation type support is "local_file_csv" or "hosted_tiledb_array')
        else:
            if self.user_annotations__type == "local_file_csv":
                dirname = self.user_annotations__local_file_csv__directory
                filename = self.user_annotations__local_file_csv__file
                if filename is not None:
                    context["messsagefn"]("Warning: --annotations-file ignored as annotations are disabled.")
                if dirname is not None:
                    context["messagefn"]("Warning: --annotations-dir ignored as annotations are disabled.")

            if self.user_annotations__ontology__enable:
                context["messagefn"](
                    "Warning: --experimental-annotations-ontology" " ignored as annotations are disabled."
                )
            if self.user_annotations__ontology__obo_location is not None:
                context["messagefn"](
                    "Warning: --experimental-annotations-ontology-obo" " ignored as annotations are disabled."
                )

    def handle_embeddings(self, context):
        self.check_attr("embeddings__names", list)
        self.check_attr("embeddings__enable_reembedding", bool)

        server_config = self.app_config.server_config
        if server_config.single_dataset__datapath:
            if self.embeddings__enable_reembedding:
                matrix_data_loader = MatrixDataLoader(
                    server_config.single_dataset__datapath, app_config=self.app_config
                )
                if matrix_data_loader.matrix_data_type != MatrixDataType.H5AD:
                    raise ConfigurationError("'enable-reembedding is only supported with H5AD files.")
                if server_config.adaptor__anndata_adaptor__backed:
                    raise ConfigurationError("enable-reembedding is not supported when run in --backed mode.")

    def handle_diffexp(self, context):
        self.check_attr("diffexp__enable", bool)
        self.check_attr("diffexp__lfc_cutoff", float)
        self.check_attr("diffexp__top_n", int)

        server_config = self.app_config.server_config
        if server_config.single_dataset__datapath:
            with server_config.matrix_data_cache_manager.data_adaptor(
                    self.tag, server_config.single_dataset__datapath, self.app_config
            ) as data_adaptor:
                if self.diffexp__enable and data_adaptor.parameters.get("diffexp_may_be_slow", False):
                    context["messagefn"](
                        "CAUTION: due to the size of your dataset, "
                        "running differential expression may take longer or fail."
                    )
