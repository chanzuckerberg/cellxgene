"""cellxgene AWS elastic beanstalk application"""

import sys
import os
import logging

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
    from backend.czi_hosted.common.config.app_config import AppConfig
    from backend.czi_hosted.app.app import Server, WSGIServer
    from backend.common.utils.data_locator import DataLocator, discover_s3_region_name
except Exception:
    logging.critical("Exception importing server modules", exc_info=True)
    sys.exit(1)

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

    # overwrite configuration for the eb app
    app_config.update_default_dataset_config(embeddings__enable_reembedding=False,)
    app_config.update_server_config(multi_dataset__allowed_matrix_types=["cxg"],)

    # complete config
    app_config.complete_config(logging.info)

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
        application.run(host=app_config.server_config.app__host, debug=debug, threaded=not debug, use_debugger=False)
    except Exception:
        logging.critical("Caught exception during initialization", exc_info=True)
        sys.exit(1)
