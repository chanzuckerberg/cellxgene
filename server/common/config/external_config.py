import os

from server.common.config.base_config import BaseConfig
from server.common.errors import ConfigurationError
from server.common.config import get_secret_key
from server.common.errors import SecretKeyRetrievalError


class ExternalConfig(BaseConfig):
    """Manages the config attribute associated with external configuration sources, such as
    environment variables or the AWS Secrets Manager."""

    def __init__(self, app_config, default_config):
        super().__init__(app_config, default_config)
        try:
            self.environment = default_config["environment"]
            self.aws_secrets_manager__region = default_config["aws_secrets_manager"]["region"]
            self.aws_secrets_manager__secrets = default_config["aws_secrets_manager"]["secrets"]

        except KeyError as e:
            raise ConfigurationError(f"Unexpected config: {str(e)}")

    @staticmethod
    def convert_value(value: str):
        """convert the external value (a string), to the most appropriate type"""
        if value == "true":
            return True
        if value == "false":
            return False
        if value == "null":
            return None
        try:
            return eval(value)
        except (NameError, SyntaxError):
            pass

        return value

    def complete_config(self, context):
        self.handle_environment(context)
        self.handle_aws_secrets_manager(context)

    def handle_environment(self, context):
        self.validate_correct_type_of_configuration_attribute("environment", list)
        for envdict in self.environment:
            name = envdict.get("name")
            if name is None:
                raise ConfigurationError("environment 'name' is missing")
            required = envdict.get("required", False)
            if type(required) != bool:
                raise ConfigurationError("environment 'required' must be a bool")
            path = envdict.get("path")
            if path is None:
                raise ConfigurationError("environment 'path' is missing")

            value = os.environ.get(name)
            if value is None:
                if required:
                    raise ConfigurationError(f"required environment variable '{name}' not set")
            else:
                value = self.convert_value(value)
                self.app_config.update_single_config_from_path_and_value(path, value)

    def handle_aws_secrets_manager(self, context):
        self.validate_correct_type_of_configuration_attribute("aws_secrets_manager__region", (type(None), str))
        self.validate_correct_type_of_configuration_attribute("aws_secrets_manager__secrets", list)

        if len(self.aws_secrets_manager__secrets) == 0:
            return

        self.validate_correct_type_of_configuration_attribute("aws_secrets_manager__region", str)

        for secret in self.aws_secrets_manager__secrets:
            secret_name = secret.get("name")
            if secret_name is None:
                raise ConfigurationError("aws_secrets_manager 'name' is missing")

            try:
                secret_dict = get_secret_key(self.aws_secrets_manager__region, secret_name)
            except SecretKeyRetrievalError as e:
                raise ConfigurationError(f"Unable to retrieve secret {secret_name}: {str(e)}")

            values = secret.get("values")
            if values is None:
                raise ConfigurationError("aws_secrets_manager 'values' is missing")
            for value in values:
                key = value.get("key")
                if key is None:
                    raise ConfigurationError(f"missing 'key' in secret values: {secret_name}")
                path = value.get("path")
                if path is None:
                    raise ConfigurationError(f"missing 'path' in secret values: {secret_name}")
                required = value.get("required", False)
                if type(required) != bool:
                    raise ConfigurationError("aws_secrets_manager 'required' must be a bool")

                secret_value = secret_dict.get(key)
                if secret_value is None:
                    if required:
                        raise ConfigurationError(f"required secret '{secret_name}:{key}' not set")
                else:
                    secret_value = self.convert_value(secret_value)
                    self.app_config.update_single_config_from_path_and_value(path, secret_value)
