import os
import sys
import warnings
from os.path import basename
from urllib.parse import urlparse

from backend.server.common.config.base_config import BaseConfig
from backend.server.common.config import DEFAULT_SERVER_PORT, BIG_FILE_SIZE_THRESHOLD
from backend.common.utils.data_locator import discover_s3_region_name
from backend.common.errors import ConfigurationError, DatasetAccessError
from backend.common.utils.utils import is_port_available, find_available_port, custom_format_warning
from backend.server.data_common.matrix_loader import MatrixDataLoader


class ServerConfig(BaseConfig):
    """Manages the config attribute associated with the server."""

    def __init__(self, app_config, default_config):
        super().__init__(app_config, default_config)

        try:
            self.app__verbose = default_config["app"]["verbose"]
            self.app__debug = default_config["app"]["debug"]
            self.app__host = default_config["app"]["host"]
            self.app__port = default_config["app"]["port"]
            self.app__open_browser = default_config["app"]["open_browser"]
            self.app__force_https = default_config["app"]["force_https"]
            self.app__flask_secret_key = default_config["app"]["flask_secret_key"]
            self.app__generate_cache_control_headers = default_config["app"]["generate_cache_control_headers"]

            self.single_dataset__datapath = default_config["single_dataset"]["datapath"]
            self.single_dataset__obs_names = default_config["single_dataset"]["obs_names"]
            self.single_dataset__var_names = default_config["single_dataset"]["var_names"]
            self.single_dataset__about = default_config["single_dataset"]["about"]
            self.single_dataset__title = default_config["single_dataset"]["title"]

            self.data_locator__s3__region_name = default_config["data_locator"]["s3"]["region_name"]

            self.adaptor__anndata_adaptor__backed = default_config["adaptor"]["anndata_adaptor"]["backed"]

            self.limits__diffexp_cellcount_max = default_config["limits"]["diffexp_cellcount_max"]
            self.limits__column_request_max = default_config["limits"]["column_request_max"]

        except KeyError as e:
            raise ConfigurationError(f"Unexpected config: {str(e)}")

        self.data_adaptor = None

    def complete_config(self, context):
        self.handle_app(context)
        self.handle_data_source()
        self.handle_data_locator()
        self.handle_adaptor()  # may depend on data_locator
        self.handle_single_dataset(context)  # may depend on adaptor
        self.handle_limits()

        self.check_config()

    def handle_app(self, context):
        self.validate_correct_type_of_configuration_attribute("app__verbose", bool)
        self.validate_correct_type_of_configuration_attribute("app__debug", bool)
        self.validate_correct_type_of_configuration_attribute("app__host", str)
        self.validate_correct_type_of_configuration_attribute("app__port", (type(None), int))
        self.validate_correct_type_of_configuration_attribute("app__open_browser", bool)
        self.validate_correct_type_of_configuration_attribute("app__force_https", bool)
        self.validate_correct_type_of_configuration_attribute("app__flask_secret_key", str)
        self.validate_correct_type_of_configuration_attribute("app__generate_cache_control_headers", bool)

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

    def handle_data_locator(self):
        self.validate_correct_type_of_configuration_attribute("data_locator__s3__region_name", (type(None), bool, str))
        if self.data_locator__s3__region_name is True:
            path = self.single_dataset__datapath

            if path.startswith("s3://"):
                region_name = discover_s3_region_name(path)
                if region_name is None:
                    raise ConfigurationError(f"Unable to discover s3 region name from {path}")
            else:
                region_name = None
            self.data_locator__s3__region_name = region_name

    def handle_data_source(self):
        self.validate_correct_type_of_configuration_attribute("single_dataset__datapath", str)

    def handle_single_dataset(self, context):
        self.validate_correct_type_of_configuration_attribute("single_dataset__datapath", (str, type(None)))
        self.validate_correct_type_of_configuration_attribute("single_dataset__title", (str, type(None)))
        self.validate_correct_type_of_configuration_attribute("single_dataset__about", (str, type(None)))
        self.validate_correct_type_of_configuration_attribute("single_dataset__obs_names", (str, type(None)))
        self.validate_correct_type_of_configuration_attribute("single_dataset__var_names", (str, type(None)))

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

    def handle_adaptor(self):
        self.validate_correct_type_of_configuration_attribute("adaptor__anndata_adaptor__backed", bool)

    def handle_limits(self):
        self.validate_correct_type_of_configuration_attribute("limits__diffexp_cellcount_max", (type(None), int))
        self.validate_correct_type_of_configuration_attribute("limits__column_request_max", (type(None), int))

    def exceeds_limit(self, limit_name, value):
        limit_value = getattr(self, "limits__" + limit_name, None)
        if limit_value is None:  # disabled
            return False
        return value > limit_value
