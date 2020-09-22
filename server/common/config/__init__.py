import logging
import os
import sys

from server.common.aws_secret_utils import get_secret_key
from server.common.data_locator import discover_s3_region_name

DEFAULT_SERVER_PORT = 5005
BIG_FILE_SIZE_THRESHOLD = 100 * 2 ** 20  # 100MB


def handle_config_from_secret(app_config):
    """Update configuration from the secret manager"""
    secret_name = os.getenv("CXG_AWS_SECRET_NAME")
    if not secret_name:
        return

    # need to find the secret manager region.
    #  1. from CXG_AWS_SECRET_REGION_NAME
    #  2. discover from dataroot location (if on s3)
    #  3. discover from config file location (if on s3)
    secret_region_name = os.getenv("CXG_AWS_SECRET_REGION_NAME")
    if secret_region_name is None:
        secret_region_name = discover_s3_region_name(app_config.multi_dataset__dataroot)
        if not secret_region_name:
            from server.eb.app import config_file
            secret_region_name = discover_s3_region_name(config_file)
    if not secret_region_name:
        logging.error("Could not determine the AWS Secret Manager region")
        sys.exit(1)

    secrets = get_secret_key(secret_region_name, secret_name)

    if not secrets:
        return

    server_attrs = (
        ("flask_secret_key", "app__flask_secret_key"),
        ("oauth_client_secret", "authentication__params_oauth__client_secret"),
    )
    default_dataset_attrs = (
        ("db_uri", "user_annotations__hosted_tiledb_array__db_uri"),
    )

    # update server configuration attributes
    for key, attr in server_attrs:
        cur_val = getattr(app_config.server_config, attr)
        if cur_val:
            continue

        # replace the attr with the secret if it is not set
        val = secrets.get(key)
        if val:
            logging.info(f"set {attr} from secret")
            app_config.update_server_config(**{attr: val})

    # update default dataset configuration attributes
    for key, attr in default_dataset_attrs:
        cur_val = getattr(app_config.default_dataset_config, attr)
        if cur_val:
            continue

        # replace the attr with the secret if it is not set
        val = secrets.get(key)
        if val:
            logging.info(f"set {attr} from secret")
            app_config.update_default_dataset_config(**{attr: val})
