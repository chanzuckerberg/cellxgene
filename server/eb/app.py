"""cellxgene AWS elastic beanstalk application"""

import sys
import os
import hashlib
import base64
from urllib.parse import urlparse
from flask import json
import logging
from flask_talisman import Talisman
from flask_cors import CORS
from server.common.aws_secret_utils import handle_config_from_secret
from server.common.errors import SecretKeyRetrievalError


if os.path.isdir("/opt/python/log"):
    # This is the standard location where Amazon EC2 instances store the application logs.
    logging.basicConfig(
        filename="/opt/python/log/app.log",
        level=logging.INFO,
        format="%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

SERVERDIR = os.path.dirname(os.path.realpath(__file__))
sys.path.append(SERVERDIR)

try:
    from server.common.app_config import AppConfig
    from server.app.app import Server
    from server.common.data_locator import DataLocator, discover_s3_region_name
except Exception:
    logging.critical("Exception importing server modules", exc_info=True)
    sys.exit(1)


class WSGIServer(Server):
    def __init__(self, app_config):
        super().__init__(app_config)

    @staticmethod
    def _before_adding_routes(app, app_config):
        script_hashes = WSGIServer.get_csp_hashes(app, app_config)
        server_config = app_config.server_config

        # add the api_base_url to the connect_src csp header.
        extra_connect_src = []
        api_base_url = server_config.get_api_base_url()
        if api_base_url:
            parse_api_base_url = urlparse(api_base_url)
            extra_connect_src = [f"{parse_api_base_url.scheme}://{parse_api_base_url.netloc}"]

        print("EXTRA CONNECT SRC", extra_connect_src)

        # This hash should be in sync with the script within
        # `client/configuration/webpack/obsoleteHTMLTemplate.html`

        # It is _very_ difficult to generate the correct hash manually,
        # consider forcing CSP to fail on the local server by intercepting the response via Requestly
        # this should print the failing script's hash to console.
        # See more here: https://github.com/chanzuckerberg/cellxgene/pull/1745
        obsolete_browser_script_hash = ["'sha256-/rmgOi/skq9MpiZxPv6lPb1PNSN+Uf4NaUHO/IjyfwM='"]
        csp = {
            "default-src": ["'self'"],
            "connect-src": ["'self'"] + extra_connect_src,
            "script-src": ["'self'", "'unsafe-eval'"]
            + obsolete_browser_script_hash + script_hashes,
            "style-src": ["'self'", "'unsafe-inline'"],
            "img-src": ["'self'", "https://cellxgene.cziscience.com", "data:"],
            "object-src": ["'none'"],
            "base-uri": ["'none'"],
            "frame-ancestors": ["'none'"],
        }

        if not app.debug:
            csp["upgrade-insecure-requests"] = ""

        if server_config.app__csp_directives:
            for k, v in server_config.app__csp_directives.items():
                if not isinstance(v, list):
                    v = [v]
                csp[k] = csp.get(k, []) + v

        # Add the web_base_url to the CORS header
        web_base_url = server_config.get_web_base_url()
        if web_base_url:
            web_base_url_parse = urlparse(web_base_url)
            allowed_origin = f"{web_base_url_parse.scheme}://{web_base_url_parse.netloc}"
            CORS(app, supports_credentials=True, origins=allowed_origin)

        Talisman(
            app, force_https=server_config.app__force_https, frame_options="DENY", content_security_policy=csp,
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
        if len(script_hashes) == 0:
            logging.error("Content security policy hashes are missing, falling back to unsafe-inline policy")

        return (script_hashes)

    @staticmethod
    def compute_inline_csp_hashes(app, app_config):
        dataset_configs = [app_config.default_dataset_config] + list(app_config.dataroot_config.values())
        hashes = []
        for dataset_config in dataset_configs:
            inline_scripts = dataset_config.app__inline_scripts
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
        script_hashes = WSGIServer.load_static_csp_hashes(app)
        script_hashes += WSGIServer.compute_inline_csp_hashes(app, app_config)
        return script_hashes


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
        app_config.update_server_config(multi_dataset__dataroot=dataroot)

    # update from secret manager
    try:
        handle_config_from_secret(app_config)
    except SecretKeyRetrievalError:
        sys.exit(1)

    # features are unsupported in the current hosted server
    app_config.update_default_dataset_config(
        embeddings__enable_reembedding=False,
    )
    app_config.update_server_config(multi_dataset__allowed_matrix_types=["cxg"],)
    app_config.complete_config(logging.info)

    if not app_config.server_config.app__flask_secret_key:
        logging.critical(
            "flask_secret_key is not provided.  Either set in config file, CXG_SECRET_KEY environment variable, "
            "or in AWS Secret Manager"
        )
        sys.exit(1)

    server = WSGIServer(app_config)

    debug = False
    application = server.app

except Exception:
    logging.critical("Caught exception during initialization", exc_info=True)
    sys.exit(1)

if app_config.is_multi_dataset():
    logging.info(f"starting server with multi_dataset__dataroot={app_config.server_config.multi_dataset__dataroot}")
else:
    logging.info(f"starting server with single_dataset__datapath={app_config.server_config.single_dataset__datapath}")

if __name__ == "__main__":
    try:
        application.run(debug=debug, threaded=not debug, use_debugger=False)
    except Exception:
        logging.critical("Caught exception during initialization", exc_info=True)
        sys.exit(1)
