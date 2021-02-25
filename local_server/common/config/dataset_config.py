import os
from os.path import splitext, isdir

from local_server.common.annotations.local_file_csv import AnnotationsLocalFile
from local_server.common.config.base_config import BaseConfig
from local_server.common.errors import ConfigurationError, OntologyLoadFailure, AnnotationsError
from local_server.compute.scanpy import get_scanpy_module
from local_server.data_common.matrix_loader import MatrixDataLoader


class DatasetConfig(BaseConfig):
    """Manages the config attribute associated with a dataset."""

    def __init__(self, tag, app_config, default_config):
        super().__init__(app_config, default_config)
        self.tag = tag
        try:
            self.app__scripts = default_config["app"]["scripts"]
            self.app__inline_scripts = default_config["app"]["inline_scripts"]
            self.app__authentication_enable = default_config["app"]["authentication_enable"]

            self.presentation__max_categories = default_config["presentation"]["max_categories"]
            self.presentation__custom_colors = default_config["presentation"]["custom_colors"]

            self.user_annotations__enable = default_config["user_annotations"]["enable"]
            self.user_annotations__type = default_config["user_annotations"]["type"]
            self.user_annotations__local_file_csv__directory = default_config["user_annotations"]["local_file_csv"][
                "directory"
            ]
            self.user_annotations__local_file_csv__file = default_config["user_annotations"]["local_file_csv"]["file"]
            self.user_annotations__ontology__enable = default_config["user_annotations"]["ontology"]["enable"]
            self.user_annotations__ontology__obo_location = default_config["user_annotations"]["ontology"][
                "obo_location"
            ]
            self.user_annotations__genesets__readonly = default_config["user_annotations"]["genesets"]["readonly"]
            self.user_annotations__local_file_csv__genesets_file = default_config["user_annotations"]["local_file_csv"][
                "genesets_file"
            ]

            self.embeddings__names = default_config["embeddings"]["names"]
            self.embeddings__enable_reembedding = default_config["embeddings"]["enable_reembedding"]

            self.diffexp__enable = default_config["diffexp"]["enable"]
            self.diffexp__lfc_cutoff = default_config["diffexp"]["lfc_cutoff"]
            self.diffexp__top_n = default_config["diffexp"]["top_n"]

        except KeyError as e:
            raise ConfigurationError(f"Unexpected config: {str(e)}")

        # The annotation object is created during complete_config and stored here.
        self.user_annotations = None

    def complete_config(self, context):
        self.handle_app()
        self.handle_presentation()
        self.handle_user_annotations(context)
        self.handle_embeddings()
        self.handle_diffexp(context)

    def get_data_adaptor(self):
        server_config = self.app_config.server_config
        if not server_config.data_adaptor:
            matrix_data_loader = MatrixDataLoader(server_config.single_dataset__datapath, app_config=self.app_config)
            server_config.data_adaptor = matrix_data_loader.open(self.app_config)

        return server_config.data_adaptor

    def handle_app(self):
        self.validate_correct_type_of_configuration_attribute("app__scripts", list)
        self.validate_correct_type_of_configuration_attribute("app__inline_scripts", list)
        self.validate_correct_type_of_configuration_attribute("app__authentication_enable", bool)

        # scripts can be string (filename) or dict (attributes). Convert string to dict.
        scripts = []
        for script in self.app__scripts:
            try:
                if isinstance(script, str):
                    scripts.append({"src": script})
                elif isinstance(script, dict) and isinstance(script["src"], str):
                    scripts.append(script)
                else:
                    raise Exception
            except Exception as e:
                raise ConfigurationError(f"Scripts must be string or a dict containing an src key: {e}")

        self.app__scripts = scripts

    def handle_presentation(self):
        self.validate_correct_type_of_configuration_attribute("presentation__max_categories", int)
        self.validate_correct_type_of_configuration_attribute("presentation__custom_colors", bool)

    def handle_user_annotations(self, context):
        self.validate_correct_type_of_configuration_attribute("user_annotations__enable", bool)
        self.validate_correct_type_of_configuration_attribute("user_annotations__type", str)
        self.validate_correct_type_of_configuration_attribute(
            "user_annotations__local_file_csv__directory", (type(None), str)
        )
        self.validate_correct_type_of_configuration_attribute(
            "user_annotations__local_file_csv__file", (type(None), str)
        )
        self.validate_correct_type_of_configuration_attribute(
            "user_annotations__local_file_csv__genesets_file", (type(None), str)
        )
        self.validate_correct_type_of_configuration_attribute("user_annotations__ontology__enable", bool)
        self.validate_correct_type_of_configuration_attribute(
            "user_annotations__ontology__obo_location", (type(None), str)
        )
        self.validate_correct_type_of_configuration_attribute("user_annotations__genesets__readonly", bool)

        if self.user_annotations__enable or not self.user_annotations__genesets__readonly:
            server_config = self.app_config.server_config
            if not self.app__authentication_enable:
                raise ConfigurationError("user annotations requires authentication to be enabled")
            if not server_config.auth.is_valid_authentication_type():
                auth_type = server_config.authentication__type
                raise ConfigurationError(f"authentication method {auth_type} is not compatible with user annotations")

        # Must always have an annotations instance to support genesets. User annotation (cell labels) are optional
        # as are writable gene sets
        if self.user_annotations__type == "local_file_csv":
            self.handle_local_file_csv_annotations(context)
        else:
            raise ConfigurationError('The only annotation type support is "local_file_csv"')

        if self.user_annotations__enable:
            if self.user_annotations__ontology__enable or self.user_annotations__ontology__obo_location:
                try:
                    self.user_annotations.load_ontology(self.user_annotations__ontology__obo_location)
                except OntologyLoadFailure as e:
                    raise ConfigurationError("Unable to load ontology terms\n" + str(e))

        self.check_annotation_config_vars_not_set(context)

    def handle_local_file_csv_annotations(self, context):
        dirname = self.user_annotations__local_file_csv__directory
        filename = self.user_annotations__local_file_csv__file
        genesets_filename = self.user_annotations__local_file_csv__genesets_file

        if dirname is not None and (filename is not None or genesets_filename is not None):
            raise ConfigurationError(
                "'user-generated-data-dir' may not be used with annotations-file' or 'genesets-file'."
            )

        if filename is not None:
            lf_name, lf_ext = splitext(filename)
            if lf_ext and lf_ext != ".csv":
                raise ConfigurationError(f"annotation file type must be .csv: {filename}")

        if genesets_filename is not None:
            lf_name, lf_ext = splitext(genesets_filename)
            if lf_ext and lf_ext != ".csv":
                raise ConfigurationError(f"genesets file type must be .csv: {genesets_filename}")

        if dirname is not None and not isdir(dirname):
            try:
                os.mkdir(dirname)
            except OSError:
                raise ConfigurationError("Unable to create directory specified by --user-generated-data-dir")

        anno_config = {
            "user-annotations": self.user_annotations__enable,
            "genesets-save": not self.user_annotations__genesets__readonly,
        }
        self.user_annotations = AnnotationsLocalFile(anno_config, dirname, filename, genesets_filename)

        # if the user has specified a fixed label file, go ahead and validate it
        # so that we can remove errors early in the process.
        server_config = self.app_config.server_config
        if server_config.single_dataset__datapath:
            data_adaptor = self.get_data_adaptor()
            if self.user_annotations__local_file_csv__file:
                data_adaptor.check_new_labels(self.user_annotations.read_labels(data_adaptor))
            if self.user_annotations__local_file_csv__genesets_file:
                try:
                    data_adaptor.check_new_genesets(self.user_annotations.read_genesets(data_adaptor, context), context)
                except (ValueError, AnnotationsError, KeyError) as e:
                    raise ConfigurationError(f"Unable to read genesets CSV file: {str(e)}") from e

    def check_annotation_config_vars_not_set(self, context):
        if self.user_annotations__type is not None:
            dirname = self.user_annotations__local_file_csv__directory
            filename = self.user_annotations__local_file_csv__file
            if not self.user_annotations__enable:
                if filename is not None:
                    context["messagefn"]("Warning: --annotations-file ignored as annotations are disabled.")
                if self.user_annotations__ontology__enable:
                    context["messagefn"](
                        "Warning: --experimental-annotations-ontology ignored as annotations are disabled."
                    )
                if self.user_annotations__ontology__obo_location is not None:
                    context["messagefn"](
                        "Warning: --experimental-annotations-ontology-obo ignored as annotations are disabled."
                    )
                if dirname is not None:
                    context["messagefn"]("Warning: --user-generated-data-dir ignored as annotations are disabled.")

    def handle_embeddings(self):
        self.validate_correct_type_of_configuration_attribute("embeddings__names", list)
        self.validate_correct_type_of_configuration_attribute("embeddings__enable_reembedding", bool)

        server_config = self.app_config.server_config
        if self.embeddings__enable_reembedding:
            if server_config.single_dataset__datapath:
                if server_config.adaptor__anndata_adaptor__backed:
                    raise ConfigurationError("enable-reembedding is not supported when run in --backed mode.")

            try:
                get_scanpy_module()
            except NotImplementedError:
                # Todo add scanpy to requirements.txt and remove this check once re-embeddings is fully supported
                raise ConfigurationError("Please install scanpy to enable UMAP re-embedding")

    def handle_diffexp(self, context):
        self.validate_correct_type_of_configuration_attribute("diffexp__enable", bool)
        self.validate_correct_type_of_configuration_attribute("diffexp__lfc_cutoff", float)
        self.validate_correct_type_of_configuration_attribute("diffexp__top_n", int)

        data_adaptor = self.get_data_adaptor()
        if self.diffexp__enable and data_adaptor.parameters.get("diffexp_may_be_slow", False):
            context["messagefn"](
                "CAUTION: due to the size of your dataset, " "running differential expression may take longer or fail."
            )
