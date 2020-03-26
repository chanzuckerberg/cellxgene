"""cellxgene AWS elastic beanstalk application"""

import sys
import os
import logging
from flask_talisman import Talisman

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
    from server.data_common.matrix_loader import MatrixDataCacheManager
    from server.common.data_locator import DataLocator
except Exception:
    logging.critical("Exception importing server modules", exc_info=True)
    sys.exit(1)


class WSGIServer(Server):
    def __init__(self, matrix_data_cache_manager, annotations, app_config):
        super().__init__(matrix_data_cache_manager, annotations, app_config)

    def _before_adding_routes(self, app_config):
        csp = {"default-src": "'self' 'unsafe-inline' 'unsafe-eval'", "img-src": ["'self'", "data:"]}
        Talisman(self.app, force_https=app_config.server__force_https, content_security_policy=csp)


try:
    app_config = AppConfig()

    dataroot = os.getenv("CXG_DATAROOT")
    config_file = os.getenv("CXG_CONFIG_FILE")

    if config_file:
        config_location = DataLocator(config_file)
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

    # features are unsupported in the current hosted server
    app_config.update(
        diffexp__enable=False,
        user_annotations__enable=False,
        embeddings__enable_reembedding=False,
        multi_dataset__allowed_matrix_types=["cxg"],
    )

    matrix_data_cache_manager = MatrixDataCacheManager()
    app_config.complete_config(matrix_data_cache_manager, logging.info)

    if not app_config.server__flask_secret_key:
        logging.critical(
            f"flask_secret_key is not provided.  Either set in config file, or in CXG_SECRET_KEY environment variable"
        )
        sys.exit(1)

    user_annotations = app_config.user_annotations

    server = WSGIServer(matrix_data_cache_manager, user_annotations, app_config)

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
