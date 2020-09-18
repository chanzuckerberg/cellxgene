import os
from os.path import splitext, isdir

import server.compute
from server.common.annotations.hosted_tiledb import AnnotationsHostedTileDB
from server.common.annotations.local_file_csv import AnnotationsLocalFile
from server.common.config.base_config import BaseConfig
from server.common.errors import ConfigurationError
from server.data_common.matrix_loader import MatrixDataLoader, MatrixDataType
from server.db.db_utils import DbUtils



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
                dc["user_annotations"][ "hosted_tiledb_array" ][ "hosted_file_directory" ]  # noqa E501

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
        if self.embeddings__enable_reembedding:
            if server_config.single_dataset__datapath:
                matrix_data_loader = MatrixDataLoader(
                    server_config.single_dataset__datapath, app_config=self.app_config
                )
                if matrix_data_loader.matrix_data_type != MatrixDataType.H5AD:
                    raise ConfigurationError("enable-reembedding is only supported with H5AD files.")
                if server_config.adaptor__anndata_adaptor__backed:
                    raise ConfigurationError("enable-reembedding is not supported when run in --backed mode.")

            try:
                server.compute.scanpy.get_scanpy_module()
            except NotImplementedError:
                raise ConfigurationError("Please install scanpy to enable UMAP re-embedding")

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
