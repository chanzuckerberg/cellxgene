"""cellxgene AWS elastic beanstalk application"""

import sys
import os
import hashlib
import base64
from flask import json
import logging
from flask_talisman import Talisman
import boto3


if os.path.isdir("/opt/python/log"):
    # This is the standard location where Amazon EC2 instances store the application logs.
    logging.basicConfig(
        filename="/opt/python/log/app.log",
        level=logging.INFO,
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
    client = session.client(service_name="secretsmanager", region_name=region_name)

    try:
        get_secret_value_response = client.get_secret_value(SecretId=secret_name)
        if "SecretString" in get_secret_value_response:
            var = get_secret_value_response["SecretString"]
            secret = json.loads(var)
            return secret.get("flask_secret_key")
    except Exception:
        logging.critical("Caught exception during get_secret_key", exc_info=True)
        sys.exit(1)

    return None


class WSGIServer(Server):
    def __init__(self, app_config):
        super().__init__(app_config)

    @staticmethod
    def _before_adding_routes(app, app_config):
        script_hashes, style_hashes = WSGIServer.get_csp_hashes(app, app_config)
        csp = {
            "default-src": ["'self'"],
            "connect-src": ["'self'"],
            "script-src": ["'self'", "'unsafe-eval'", "'unsafe-inline'"] + script_hashes,
            "style-src": ["'self'", "'unsafe-inline'"] + style_hashes,
            "img-src": ["'self'", "data:"],
            "object-src": ["'none'"],
            "base-uri": ["'none'"],
            "frame-ancestors": ["'none'"],
        }

        if not app.debug:
            csp["upgrade-insecure-requests"] = ""

        if app_config.server__csp_directives:
            for k, v in app_config.server__csp_directives.items():
                if not isinstance(v, list):
                    v = [v]
                csp[k] = csp.get(k, []) + v

        Talisman(
            app, force_https=app_config.server__force_https, frame_options="DENY", content_security_policy=csp,
        )

    @staticmethod
    def load_static_csp_hashes(app):
        csp_hashes = None
        try:
            with app.open_resource("../common/web/csp-hashes.json") as f:
                csp_hashes = json.load(f)
        except FileNotFoundError:
            pass
        if not isinstance(csp_hashes, dict):
            csp_hashes = {}
        script_hashes = [f"'{hash}'" for hash in csp_hashes.get("script-hashes", [])]
        style_hashes = [f"'{hash}'" for hash in csp_hashes.get("style-hashes", [])]

        if len(script_hashes) == 0 or len(style_hashes) == 0:
            logging.error("Content security policy hashes are missing, falling back to unsafe-inline policy")

        return (script_hashes, style_hashes)

    @staticmethod
    def compute_inline_scp_hashes(app, app_config):
        inline_scripts = app_config.server__inline_scripts
        hashes = []
        for script in inline_scripts:
            with app.open_resource(f"../common/web/templates/{script}") as f:
                content = f.read()
                # we use jinja2 template include, which trims final newline if present.
                if content[-1] == 0x0A:
                    content = content[0:-1]
                hash = base64.b64encode(hashlib.sha256(content).digest())
                hashes.append(f"'sha256-{hash.decode('utf-8')}'")
        return hashes

    @staticmethod
    def get_csp_hashes(app, app_config):
        script_hashes, style_hashes = WSGIServer.load_static_csp_hashes(app)
        script_hashes += WSGIServer.compute_inline_scp_hashes(app, app_config)
        return (script_hashes, style_hashes)


try:
    app_config = AppConfig()

    has_config = False
    # config file: look first for "config.yaml" in the current working directory
    config_file = "config.yaml"
    config_location = DataLocator(config_file)
    if config_location.exists():
        with config_location.local_handle() as lh:
            logging.info(f"Configuration from {config_file}")
            app_config.update_from_config_file(lh)
            has_config = True

    else:
        # config file: second, use the CXG_CONFIG_FILE
        config_file = os.getenv("CXG_CONFIG_FILE")
        if config_file:
            region_name = discover_s3_region_name(config_file)
            config_location = DataLocator(config_file, region_name)
            if config_location.exists():
                with config_location.local_handle() as lh:
                    logging.info(f"Configuration from {config_file}")
                    app_config.update_from_config_file(lh)
                    has_config = True
            else:
                logging.critical(f"Configuration file not found {config_file}")
                sys.exit(1)

    if not has_config:
        logging.critical("No config file found")
        sys.exit(1)

    dataroot = os.getenv("CXG_DATAROOT")
    if dataroot:
        logging.info("Configuration from CXG_DATAROOT")
        app_config.update(multi_dataset__dataroot=dataroot)

    secret_name = os.getenv("CXG_AWS_SECRET_NAME")
    if secret_name:
        # need to find the secret manager region.
        #  1. from CXG_AWS_SECRET_REGION_NAME
        #  2. discover from dataroot location (if on s3)
        #  3. discover from config file location (if on s3)
        secret_region_name = os.getenv("CXG_AWS_SECRET_REGION_NAME")
        if secret_region_name is None:
            secret_region_name = discover_s3_region_name(app_config.multi_dataset__dataroot)
            if not secret_region_name:
                secret_region_name = discover_s3_region_name(config_file)
        if not secret_region_name:
            logging.error("Could not determine the AWS Secret Manager region")
            sys.exit(1)

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
            "flask_secret_key is not provided.  Either set in config file, CXG_SECRET_KEY environment variable, "
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
