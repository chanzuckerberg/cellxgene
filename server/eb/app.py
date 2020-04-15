"""cellxgene AWS elastic beanstalk application"""

import sys
import os
import logging
from flask_talisman import Talisman
import boto3
import json

if os.path.isdir("/opt/python/log"):
    # This is the standard location where Amazon EC2 instances store the application logs.
    logging.basicConfig(
        filename="/opt/python/log/app.log",
        level=logging.DEBUG,
        format="%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

# echo the logs to stdout.  Useful for local testing
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

SERVERDIR = os.path.dirname(os.path.realpath(__file__))
sys.path.append(SERVERDIR)

try:
    from server.common.app_config import AppConfig
    from server.app.app import Server
    from server.common.data_locator import DataLocator, discover_s3_region_name
except Exception:
    logging.critical("Exception importing server modules", exc_info=True)
    sys.exit(1)


def get_flask_secret_key(region_name, secret_name):
    session = boto3.session.Session()
    client = session.client(
        service_name='secretsmanager',
        region_name=region_name
    )

    try:
        get_secret_value_response = client.get_secret_value(
            SecretId=secret_name
        )
        if 'SecretString' in get_secret_value_response:
            var = get_secret_value_response['SecretString']
            secret = json.loads(var)
            return secret.get("flask_secret_key")
    except Exception:
        logging.critical("Caught exception during get_secret_key", exc_info=True)
        sys.exit(1)

    return None


class WSGIServer(Server):
    def __init__(self, app_config):
        super().__init__(app_config)

    def _before_adding_routes(self, app_config):
        csp = {"default-src": "'self' 'unsafe-inline' 'unsafe-eval'", "img-src": ["'self'", "data:"]}
        Talisman(self.app, force_https=app_config.server__force_https, content_security_policy=csp)


try:
    app_config = AppConfig()

    dataroot = os.getenv("CXG_DATAROOT")
    config_file = os.getenv("CXG_CONFIG_FILE")

    secret_name = os.getenv("CXG_AWS_SECRET_NAME")
    secret_region_name = os.getenv("CXG_AWS_SECRET_REGION_NAME")

    if config_file:
        region_name = discover_s3_region_name(config_file)
        config_location = DataLocator(config_file, region_name)
        if config_location.exists():
            with config_location.local_handle() as lh:
                logging.info(f"Configuration from {config_file}")
                app_config.update_from_config_file(lh)
        else:
            logging.critical(f"Configuration file not found {config_file}")
            sys.exit(1)
    else:
        # no config file specified, try "config.yaml" in the current working directory
        config_file = "config.yaml"
        config_location = DataLocator(config_file)
        if config_location.exists():
            with config_location.local_handle() as lh:
                logging.info(f"Configuration from {config_file}")
                app_config.update_from_config_file(lh)

    if dataroot:
        logging.info(f"Configuration from CXG_DATAROOT")
        app_config.update(multi_dataset__dataroot=dataroot)

    if secret_name:
        if secret_region_name is None:
            secret_region_name = discover_s3_region_name(app_config.multi_dataset__dataroot)
            if not secret_region_name:
                logging.error(f"Expected to discover the s3 region name from {app_config.multi_dataset__dataroot}")
        flask_secret_key = get_flask_secret_key(secret_region_name, secret_name)
        app_config.update(server__flask_secret_key=flask_secret_key)

    # features are unsupported in the current hosted server
    app_config.update(
        user_annotations__enable=False,
        embeddings__enable_reembedding=False,
        multi_dataset__allowed_matrix_types=["cxg"],
    )

    app_config.complete_config(logging.info)

    if not app_config.server__flask_secret_key:
        logging.critical(
            f"flask_secret_key is not provided.  Either set in config file, CXG_SECRET_KEY environment variable, "
            "or in AWS Secret Manager"
        )
        sys.exit(1)

    user_annotations = app_config.user_annotations

    server = WSGIServer(app_config)

    debug = False
    application = server.app

except Exception:
    logging.critical("Caught exception during initialization", exc_info=True)
    sys.exit(1)

if app_config.multi_dataset__dataroot:
    logging.info(f"starting server with multi_dataset__dataroot={app_config.multi_dataset__dataroot}")
elif app_config.single_dataset__datapath:
    logging.info(f"starting server with single_dataset__datapath={app_config.single_dataset__datapath}")

if __name__ == "__main__":
    try:
        application.run(debug=debug, threaded=not debug, use_debugger=False)
    except Exception:
        logging.critical("Caught exception during initialization", exc_info=True)
        sys.exit(1)
