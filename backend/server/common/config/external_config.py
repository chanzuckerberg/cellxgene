import os

from backend.server.common.config.base_config import BaseConfig
from backend.common.errors import ConfigurationError
from backend.common.utils.type_conversion_utils import convert_string_to_value


class ExternalConfig(BaseConfig):
    """Manages the config attribute associated with external configuration sources, such as
    environment variables."""

    def __init__(self, app_config, default_config):
        super().__init__(app_config, default_config)
        try:
            self.environment = default_config["environment"]

        except KeyError as e:
            raise ConfigurationError(f"Unexpected config: {str(e)}")

    def complete_config(self, context):
        self.handle_environment(context)

    def handle_environment(self, context):
        """For each environment variable defined, get the value (if it is set),
        and set the specified config parameter"""
        self.validate_correct_type_of_configuration_attribute("environment", list)
        for envdict in self.environment:
            name = envdict.get("name")
            if name is None:
                raise ConfigurationError("environment: 'name' is missing")
            required = envdict.get("required", False)
            if type(required) != bool:
                raise ConfigurationError("environment: 'required' must be a bool")
            path = envdict.get("path")
            if path is None:
                raise ConfigurationError("environment: 'path' is missing")

            value = os.environ.get(name)
            if value is None:
                if required:
                    raise ConfigurationError(f"required environment variable '{name}' not set")
            else:
                value = convert_string_to_value(value)
                self.app_config.update_single_config_from_path_and_value(path, value)
