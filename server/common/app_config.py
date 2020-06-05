from server import __version__ as cellxgene_version
from flatten_dict import flatten, unflatten
import os
from os.path import splitext, basename, isdir
import sys
from urllib.parse import urlparse, quote_plus
import yaml
import copy

from server.common.default_config import get_default_config
from server.common.errors import ConfigurationError, DatasetAccessError, OntologyLoadFailure
from server.data_common.matrix_loader import MatrixDataLoader, MatrixDataCacheManager, MatrixDataType
from server.common.utils import find_available_port, is_port_available
import warnings
from server.common.annotations import AnnotationsLocalFile
from server.common.utils import custom_format_warning
import server.compute.diffexp_cxg as diffexp_tiledb
from server.common.data_locator import discover_s3_region_name

DEFAULT_SERVER_PORT = int(os.environ.get("CXG_SERVER_PORT", "5005"))
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
    def __init__(self):

        self.default_config = get_default_config()
        self.attr_checked = {k: False for k in self.__mapping(self.default_config).keys()}

        dc = self.default_config
        try:
            self.server__verbose = dc["server"]["verbose"]
            self.server__debug = dc["server"]["debug"]
            self.server__host = dc["server"]["host"]
            self.server__port = dc["server"]["port"]
            self.server__scripts = dc["server"]["scripts"]
            self.server__inline_scripts = dc["server"]["inline_scripts"]
            self.server__open_browser = dc["server"]["open_browser"]
            self.server__about_legal_tos = dc["server"]["about_legal_tos"]
            self.server__about_legal_privacy = dc["server"]["about_legal_privacy"]
            self.server__force_https = dc["server"]["force_https"]
            self.server__flask_secret_key = dc["server"]["flask_secret_key"]
            self.server__generate_cache_control_headers = dc["server"]["generate_cache_control_headers"]
            self.server__server_timing_headers = dc["server"]["server_timing_headers"]
            self.server__csp_directives = dc["server"]["csp_directives"]

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

            self.user_annotations__enable = dc["user_annotations"]["enable"]
            self.user_annotations__type = dc["user_annotations"]["type"]
            self.user_annotations__local_file_csv__directory = dc["user_annotations"]["local_file_csv"]["directory"]
            self.user_annotations__local_file_csv__file = dc["user_annotations"]["local_file_csv"]["file"]
            self.user_annotations__ontology__enable = dc["user_annotations"]["ontology"]["enable"]
            self.user_annotations__ontology__obo_location = dc["user_annotations"]["ontology"]["obo_location"]

            self.presentation__max_categories = dc["presentation"]["max_categories"]
            self.presentation__custom_colors = dc["presentation"]["custom_colors"]

            self.embeddings__names = dc["embeddings"]["names"]
            self.embeddings__enable_reembedding = dc["embeddings"]["enable_reembedding"]

            self.diffexp__enable = dc["diffexp"]["enable"]
            self.diffexp__lfc_cutoff = dc["diffexp"]["lfc_cutoff"]
            self.diffexp__top_n = dc["diffexp"]["top_n"]
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

        # The annotation object is created during complete_config and stored here.
        self.user_annotations = None

        # The matrix data cache manager is created during the complete_config and stored here.
        self.matrix_data_cache_manager = None

        # Set to true when config_completed is called
        self.is_completed = False

    def check_config(self):
        if not self.is_completed:
            raise ConfigurationError("The configuration has not been completed")
        mapping = self.__mapping(self.default_config)
        for key in mapping.keys():
            if not self.attr_checked[key]:
                raise ConfigurationError(f"The attr '{key}' has not been checked")

    def __mapping(self, config):
        """Create a mapping from attribute names to (location in the config tree, value)"""
        dc = copy.deepcopy(config)
        mapping = {}

        # special cases where the value could be a dict.
        # If its value is not None, the entry is added to the mapping, and not included
        # in the flattening below.
        dictval_cases = [
            ("adaptor", "cxg_adaptor", "tiledb_ctx"),
            ("server", "csp_directives"),
            ("multi_dataset", "dataroot"),
        ]
        for dictval_case in dictval_cases:
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

    def update_from_config_file(self, config_file):
        with open(config_file) as fyaml:
            config = yaml.load(fyaml, Loader=yaml.FullLoader)

        mapping = self.__mapping(config)
        for attr, (key, value) in mapping.items():
            if not hasattr(self, attr):
                raise ConfigurationError(f"Unknown key from config file: {key}")
            try:
                setattr(self, attr, value)
            except KeyError:
                raise ConfigurationError(f"Unable to set config attribute: {key}")

            self.attr_checked[attr] = False

        self.is_completed = False

    def write_config(self, config_file):
        """output the config to a yaml file"""
        mapping = self.__mapping(self.default_config)
        for attrname in mapping.keys():
            mapping[attrname] = getattr(self, attrname)
        config = unflatten(mapping, splitter=lambda key: key.split("__"))
        yaml.dump(config, open(config_file, "w"))

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

        self.is_completed = False

    def changes_from_default(self):
        """Return all the attribute that are different from the default"""
        mapping = self.__mapping(self.default_config)
        diff = []
        for attrname, (key, defval) in mapping.items():
            curval = getattr(self, attrname)
            if curval != defval:
                diff.append((attrname, curval, defval))
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

        self.handle_server(context)
        self.handle_adaptor(context)
        self.handle_data_locator(context)
        self.handle_adaptor(context)  # may depend on data_locator
        self.handle_presentation(context)
        self.handle_single_dataset(context)  # may depend on adaptor
        self.handle_multi_dataset(context)  # may depend on adaptor
        self.handle_user_annotations(context)
        self.handle_embeddings(context)
        self.handle_diffexp(context)
        self.handle_limits(context)

        self.is_completed = True
        self.check_config()

    def __check_attr(self, attrname, vtype):
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

    def handle_server(self, context):
        self.__check_attr("server__verbose", bool)
        self.__check_attr("server__debug", bool)
        self.__check_attr("server__host", str)
        self.__check_attr("server__port", (type(None), int))
        self.__check_attr("server__scripts", list)
        self.__check_attr("server__inline_scripts", list)
        self.__check_attr("server__open_browser", bool)
        self.__check_attr("server__force_https", bool)
        self.__check_attr("server__flask_secret_key", (type(None), str))
        self.__check_attr("server__generate_cache_control_headers", bool)
        self.__check_attr("server__about_legal_tos", (type(None), str))
        self.__check_attr("server__about_legal_privacy", (type(None), str))
        self.__check_attr("server__server_timing_headers", bool)
        self.__check_attr("server__csp_directives", (type(None), dict))

        if self.server__port:
            if not is_port_available(self.server__host, self.server__port):
                raise ConfigurationError(
                    f"The port selected {self.server__port} is in use, please configure an open port."
                )
        else:
            self.server__port = find_available_port(self.server__host, DEFAULT_SERVER_PORT)

        if self.server__debug:
            context["messagefn"]("in debug mode, setting verbose=True and open_browser=False")
            self.server__verbose = True
            self.server__open_browser = False
        else:
            warnings.formatwarning = custom_format_warning

        if not self.server__verbose:
            sys.tracebacklimit = 0

        # secret key:
        #   first, from CXG_SECRET_KEY environment variable
        #   second, from config file
        self.server__flask_secret_key = os.environ.get("CXG_SECRET_KEY", self.server__flask_secret_key)

        # CSP Directives are a dict of string: list(string) or string: string
        if self.server__csp_directives is not None:
            for k, v in self.server__csp_directives.items():
                if not isinstance(k, str):
                    raise ConfigurationError("CSP directive names must be a string.")
                if isinstance(v, list):
                    for policy in v:
                        if not isinstance(policy, str):
                            raise ConfigurationError("CSP directive value must be a string or list of strings.")
                elif not isinstance(v, str):
                    raise ConfigurationError("CSP directive value must be a string or list of strings.")

        # scripts can be string (filename) or dict (attributes).   Convert string to dict.
        scripts = []
        for s in self.server__scripts:
            if isinstance(s, str):
                scripts.append({"src": s})
            elif isinstance(s, dict) and isinstance(s["src"], str):
                scripts.append(s)
            else:
                raise ConfigurationError("Scripts must be string or dict")
        self.server__scripts = scripts

    def handle_data_locator(self, context):
        self.__check_attr("data_locator__s3__region_name", (type(None), bool, str))
        if self.data_locator__s3__region_name is True:
            path = self.single_dataset__datapath or self.multi_dataset__dataroot
            if type(path) == dict:
                # if multi_dataset__dataroot is a dict, then use the first key
                # that is in s3.   NOTE:  it is not supported to have dataroots
                # in different regions.
                paths = path.values()
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

    def handle_presentation(self, context):
        self.__check_attr("presentation__max_categories", int)
        self.__check_attr("presentation__custom_colors", bool)

    def handle_single_dataset(self, context):
        self.__check_attr("single_dataset__datapath", (str, type(None)))
        self.__check_attr("single_dataset__title", (str, type(None)))
        self.__check_attr("single_dataset__about", (str, type(None)))
        self.__check_attr("single_dataset__obs_names", (str, type(None)))
        self.__check_attr("single_dataset__var_names", (str, type(None)))

        if self.single_dataset__datapath is None:
            if self.multi_dataset__dataroot is None:
                # TODO:  change the error message once dataroot is fully supported
                raise ConfigurationError("missing datapath")
            return
        else:
            if self.multi_dataset__dataroot is not None:
                raise ConfigurationError("must supply only one of datapath or dataroot")

        # create the matrix data cache manager:
        if self.matrix_data_cache_manager is None:
            self.matrix_data_cache_manager = MatrixDataCacheManager(max_cached=1, timelimit_s=None)

        # preload this data set
        matrix_data_loader = MatrixDataLoader(self.single_dataset__datapath, app_config=self)
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
        self.__check_attr("multi_dataset__dataroot", (type(None), dict, str))
        self.__check_attr("multi_dataset__index", (type(None), bool, str))
        self.__check_attr("multi_dataset__allowed_matrix_types", list)
        self.__check_attr("multi_dataset__matrix_cache__max_datasets", int)
        self.__check_attr("multi_dataset__matrix_cache__timelimit_s", (type(None), int, float))

        if self.multi_dataset__dataroot is None:
            return

        if type(self.multi_dataset__dataroot) == str:
            self.multi_dataset__dataroot = dict(d=self.multi_dataset__dataroot)

        for key in self.multi_dataset__dataroot.keys():
            # sanity check for well formed keys
            if type(key) != str:
                raise ConfigurationError(f"error in multi_dataset__dataroot {key}")
            if quote_plus(key) != key:
                raise ConfigurationError(f"error in multi_dataset__dataroot {key}")
            if os.path.split(os.path.normpath(key))[-1] != key:
                raise ConfigurationError(f"error in multi_dataset__dataroot {key}")

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

    def handle_user_annotations(self, context):
        self.__check_attr("user_annotations__enable", bool)
        self.__check_attr("user_annotations__type", str)
        self.__check_attr("user_annotations__local_file_csv__directory", (type(None), str))
        self.__check_attr("user_annotations__local_file_csv__file", (type(None), str))
        self.__check_attr("user_annotations__ontology__enable", bool)
        self.__check_attr("user_annotations__ontology__obo_location", (type(None), str))

        if self.user_annotations__enable:
            # TODO, replace this with a factory pattern once we have more than one way
            # to do annotations.  currently only local_file_csv
            if self.user_annotations__type != "local_file_csv":
                raise ConfigurationError('The only annotation type support is "local_file_csv"')

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
            if self.single_dataset__datapath and self.user_annotations__local_file_csv__file:
                with self.matrix_data_cache_manager.data_adaptor(self.single_dataset__datapath, self) as data_adaptor:
                    data_adaptor.check_new_labels(self.user_annotations.read_labels(data_adaptor))

            if self.user_annotations__ontology__enable or self.user_annotations__ontology__obo_location:
                try:
                    self.user_annotations.load_ontology(self.user_annotations__ontology__obo_location)
                except OntologyLoadFailure as e:
                    raise ConfigurationError("Unable to load ontology terms\n" + str(e))

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
        self.__check_attr("embeddings__names", list)
        self.__check_attr("embeddings__enable_reembedding", bool)

        if self.single_dataset__datapath:
            if self.embeddings__enable_reembedding:
                matrix_data_loader = MatrixDataLoader(self.single_dataset__datapath, app_config=self)
                if matrix_data_loader.matrix_data_type() != MatrixDataType.H5AD:
                    raise ConfigurationError("'enable-reembedding is only supported with H5AD files.")
                if self.adaptor__anndata_adaptor__backed:
                    raise ConfigurationError("enable-reembedding is not supported when run in --backed mode.")

    def handle_diffexp(self, context):
        self.__check_attr("diffexp__enable", bool)
        self.__check_attr("diffexp__lfc_cutoff", float)
        self.__check_attr("diffexp__top_n", int)
        self.__check_attr("diffexp__alg_cxg__max_workers", (str, int))
        self.__check_attr("diffexp__alg_cxg__cpu_multiplier", int)
        self.__check_attr("diffexp__alg_cxg__target_workunit", int)

        if self.single_dataset__datapath:
            with self.matrix_data_cache_manager.data_adaptor(self.single_dataset__datapath, self) as data_adaptor:
                if self.diffexp__enable and data_adaptor.parameters.get("diffexp_may_be_slow", False):
                    context["messagefn"](
                        "CAUTION: due to the size of your dataset, "
                        "running differential expression may take longer or fail."
                    )

        max_workers = self.diffexp__alg_cxg__max_workers
        cpu_multiplier = self.diffexp__alg_cxg__cpu_multiplier
        cpu_count = os.cpu_count()
        max_workers = min(max_workers, cpu_multiplier * cpu_count)
        diffexp_tiledb.set_config(max_workers, self.diffexp__alg_cxg__target_workunit)

    def handle_adaptor(self, context):
        # cxg
        self.__check_attr("adaptor__cxg_adaptor__tiledb_ctx", dict)
        regionkey = "vfs.s3.region"
        if regionkey not in self.adaptor__cxg_adaptor__tiledb_ctx:
            if type(self.data_locator__s3__region_name) == str:
                self.adaptor__cxg_adaptor__tiledb_ctx[regionkey] = self.data_locator__s3__region_name

        from server.data_cxg.cxg_adaptor import CxgAdaptor

        CxgAdaptor.set_tiledb_context(self.adaptor__cxg_adaptor__tiledb_ctx)

        # anndata
        self.__check_attr("adaptor__anndata_adaptor__backed", bool)

    def handle_limits(self, context):
        self.__check_attr("limits__diffexp_cellcount_max", (type(None), int))
        self.__check_attr("limits__column_request_max", (type(None), int))

    def get_title(self, data_adaptor):
        return self.single_dataset__title if self.single_dataset__title else data_adaptor.get_title()

    def get_about(self, data_adaptor):
        return self.single_dataset__about if self.single_dataset__about else data_adaptor.get_about()

    def get_client_config(self, data_adaptor, annotation=None):
        """
        Return the configuration as required by the /config REST route
        """

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
        library_versions["cellxgene"] = cellxgene_version

        # links
        links = {"about-dataset": about}

        # parameters
        parameters = {
            "layout": self.embeddings__names,
            "max-category-items": self.presentation__max_categories,
            "obs_names": self.single_dataset__obs_names,
            "var_names": self.single_dataset__var_names,
            "diffexp_lfc_cutoff": self.diffexp__lfc_cutoff,
            "backed": self.adaptor__anndata_adaptor__backed,
            "disable-diffexp": not self.diffexp__enable,
            "enable-reembedding": self.embeddings__enable_reembedding,
            "annotations": False,
            "annotations_file": None,
            "annotations_dir": None,
            "annotations_cell_ontology_enabled": False,
            "annotations_cell_ontology_obopath": None,
            "annotations_cell_ontology_terms": None,
            "custom_colors": self.presentation__custom_colors,
            "diffexp-may-be-slow": False,
            "about_legal_tos": self.server__about_legal_tos,
            "about_legal_privacy": self.server__about_legal_privacy,
        }

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
        config["limits"] = {
            "column_request_max": self.limits__column_request_max,
            "diffexp_cellcount_max": self.limits__diffexp_cellcount_max,
        }

        return c

    def exceeds_limit(self, limit_name, value):
        limit_value = getattr(self, "limits__" + limit_name, None)
        if limit_value is None:  # disabled
            return False
        return value > limit_value
