# -*- coding: utf-8 -*-


from server import __version__ as cellxgene_version
from flatten_dict import flatten
from os import mkdir, environ
from os.path import splitext, basename, isdir
import sys
from urllib.parse import urlparse
import yaml

from server.common.default_config import get_default_config
from server.common.errors import ConfigurationError, DatasetAccessError, OntologyLoadFailure
from server.data_common.matrix_loader import MatrixDataLoader, MatrixDataCacheManager, MatrixDataType
from server.common.utils import find_available_port, is_port_available
import warnings
from server.common.annotations import AnnotationsLocalFile
from server.common.utils import custom_format_warning

DEFAULT_SERVER_PORT = int(environ.get("CXG_SERVER_PORT", "5005"))
# anything bigger than this will generate a special message
BIG_FILE_SIZE_THRESHOLD = 100 * 2 ** 20  # 100MB

""" Default limits for requests """
Default_Limits = {
    # Max number of columns that may be requested for /annotations or /data routes.
    # This is a simplistic means of preventing excess resource consumption (eg,
    # requesting the entire X matrix in one request) or other DoS style attacks/errors.
    # Set to None to disable check.
    "column_request_max": 32,
    # Max number of cells that will be accepted for differential expression.
    # Set to None to disable the check.
    "diffexp_cellcount_max": None,  # None is disabled
}


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

        dc = self.default_config
        try:
            self.server__verbose = dc["server"]["verbose"]
            self.server__debug = dc["server"]["debug"]
            self.server__host = dc["server"]["host"]
            self.server__port = dc["server"]["port"]
            self.server__scripts = dc["server"]["scripts"]
            self.server__open_browser = dc["server"]["open_browser"]
            self.server__about_legal_tos = dc["server"]["about_legal_tos"]
            self.server__about_legal_privacy = dc["server"]["about_legal_privacy"]
            self.server__force_https = dc["server"]["force_https"]
            self.server__flask_secret_key = dc["server"]["flask_secret_key"]

            self.multi_dataset__dataroot = dc["multi_dataset"]["dataroot"]
            self.multi_dataset__index = dc["multi_dataset"]["index"]
            self.multi_dataset__allowed_matrix_types = dc["multi_dataset"]["allowed_matrix_types"]
            self.multi_dataset__matrix_cache__max_datasets = dc["multi_dataset"]["matrix_cache"]["max_datasets"]

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

            self.embeddings__names = dc["embeddings"]["names"]
            self.embeddings__enable_reembedding = dc["embeddings"]["enable_reembedding"]

            self.diffexp__enable = dc["diffexp"]["enable"]
            self.diffexp__lfc_cutoff = dc["diffexp"]["lfc_cutoff"]

            self.adaptor__cxg_adaptor__tiledb_ctx = dc["adaptor"]["cxg_adaptor"]["tiledb_ctx"]
            self.adaptor__anndata_adaptor__backed = dc["adaptor"]["anndata_adaptor"]["backed"]

        except KeyError as e:
            raise ConfigurationError(f"Unexpected config: {str(e)}")

        # Used for various limits, eg, size of requests.  Not currently configurable.
        self.limits = Default_Limits

        # The annotation object is created during complete_config and stored here.
        self.user_annotations = None

        # Set to true when config_completed is called
        self.is_completed = False

    def update_from_config_file(self, config_file):
        with open(config_file) as fyaml:
            config = yaml.load(fyaml, Loader=yaml.FullLoader)

        # special case for tiledb_ctx whose value is a dict, and cannot
        # be handled by the flattening below
        if config.get("adaptor", {}).get("cxg_adaptor", {}).get("tiledb_ctx"):
            value = config["adaptor"]["cxg_adaptor"]["tiledb_ctx"]
            self.adaptor__cxg_adaptor__tiledb_ctx = value
            del config["adaptor"]["cxg_adaptor"]["tiledb_ctx"]

        flat_config = flatten(config)
        for key, value in flat_config.items():
            # name of the attribute
            attr = "__".join(key)
            if not hasattr(self, attr):
                raise ConfigurationError(f"Unknown key from config file: {key}")
            try:
                setattr(self, attr, value)
            except KeyError:
                raise ConfigurationError(f"Unable to set config attribute: {key}")

        self.is_completed = False

    def update(self, **kw):
        for key, value in kw.items():
            if not hasattr(self, key):
                raise ConfigurationError(f"unknown config parameter {key}.")
            try:
                setattr(self, key, value)
            except KeyError:
                raise ConfigurationError(f"Unable to set config parameter {key}.")

        self.is_completed = False

    def complete_config(self, matrix_data_cache_manager=None, messagefn=None):
        """The configure options are checked, and any additional setup based on the config
        parameters is done"""

        if matrix_data_cache_manager is None:
            matrix_data_cache_manager = MatrixDataCacheManager()
        if messagefn is None:

            def noop(message):
                pass

            messagefn = noop

        # TODO: to give better error messages we can add a mapping between where each config
        # attribute originated (e.g. command line argument or config file), then in the error
        # messages we can give correct context for attributes with bad value.
        context = dict(matrix_cache=matrix_data_cache_manager, messagefn=messagefn)

        self.handle_server(context)
        self.handle_single_dataset(context)
        self.handle_multi_dataset(context)
        self.handle_user_annotations(context)
        self.handle_embeddings(context)
        self.handle_diffexp(context)
        self.handle_adaptor(context)

        self.is_completed = True

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

    def handle_server(self, context):
        self.__check_attr("server__verbose", bool)
        self.__check_attr("server__debug", bool)
        self.__check_attr("server__host", str)
        self.__check_attr("server__port", (type(None), int))
        self.__check_attr("server__scripts", (list, tuple))
        self.__check_attr("server__open_browser", bool)
        self.__check_attr("server__force_https", bool)
        self.__check_attr("server__flask_secret_key", (type(None), str))

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
        self.server__flask_secret_key = environ.get("CXG_SECRET_KEY", self.server__flask_secret_key)

    def handle_presentation(self, context):
        self.__check_attr("presentation__max_categories", int)

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

        # preload this data set
        matrix_data_loader = MatrixDataLoader(self.single_dataset__datapath)
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
        self.__check_attr("multi_dataset__dataroot", (type(None), str))
        self.__check_attr("multi_dataset__index", (type(None), bool, str))
        self.__check_attr("multi_dataset__allowed_matrix_types", (tuple, list))
        self.__check_attr("multi_dataset__matrix_cache__max_datasets", int)

        if self.multi_dataset__dataroot is None:
            return

        # error checking
        for mtype in self.multi_dataset__allowed_matrix_types:
            try:
                MatrixDataType(mtype)
            except ValueError:
                raise ConfigurationError(f'Invalid matrix type in "allowed_matrix_types": {mtype}')

        # matrix cache
        MatrixDataCacheManager.set_max_datasets(self.multi_dataset__matrix_cache__max_datasets)

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
                    mkdir(dirname)
                except OSError:
                    raise ConfigurationError("Unable to create directory specified by --annotations-dir")

            self.user_annotations = AnnotationsLocalFile(dirname, filename)

            # if the user has specified a fixed label file, go ahead and validate it
            # so that we can remove errors early in the process.
            if self.single_dataset__datapath and self.user_annotations__local_file_csv__file:
                with context["matrix_cache"].data_adaptor(self.single_dataset__datapath, self) as data_adaptor:
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
        self.__check_attr("embeddings__names", (list, tuple))
        self.__check_attr("embeddings__enable_reembedding", bool)

        if self.single_dataset__datapath:
            if self.embeddings__enable_reembedding:
                matrix_data_loader = MatrixDataLoader(self.single_dataset__datapath)
                if matrix_data_loader.matrix_data_type() != MatrixDataType.H5AD:
                    raise ConfigurationError("'enable-reembedding is only supported with H5AD files.")
                if self.adaptor__anndata_adaptor__backed:
                    raise ConfigurationError("enable-reembedding is not supported when run in --backed mode.")

    def handle_diffexp(self, context):
        self.__check_attr("diffexp__enable", bool)
        self.__check_attr("diffexp__lfc_cutoff", float)

        if self.single_dataset__datapath:
            with context["matrix_cache"].data_adaptor(self.single_dataset__datapath, self) as data_adaptor:
                if self.diffexp__enable and data_adaptor.parameters.get("diffexp_may_be_slow", False):
                    context["messagefn"](
                        f"CAUTION: due to the size of your dataset, "
                        f"running differential expression may take longer or fail."
                    )

    def handle_adaptor(self, context):
        # cxg
        self.__check_attr("adaptor__cxg_adaptor__tiledb_ctx", dict)
        from server.data_cxg.cxg_adaptor import CxgAdaptor

        CxgAdaptor.set_tiledb_context(self.adaptor__cxg_adaptor__tiledb_ctx)

        # anndata
        self.__check_attr("adaptor__anndata_adaptor__backed", bool)

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

        if not self.is_completed:
            raise ConfigurationError("The configuration has not been completed")

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
        config["limits"] = self.limits

        return c

    def exceeds_limit(self, limit_name, value):
        limit_value = self.limits.get(limit_name, None)
        if limit_value is None:  # disabled
            return False
        return value > limit_value
