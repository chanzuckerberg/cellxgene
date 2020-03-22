"""cellxgene AWS elastic beanstalk application"""

import sys
import os
import logging

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
except Exception:
    logging.critical("Exception importing server modules", exc_info=True)
    sys.exit(1)

try:
    dataroot = os.getenv("CXG_DATAROOT")
    app_config = AppConfig()

    config_file = "config.yaml"
    if os.path.exists(config_file):
        logging.info(f"Configuration from {config_file}")
        app_config.update_from_config_file(config_file)

    if dataroot:
        logging.info(f"Configuration from CXG_DATAROOT")
        app_config.update(
            multi_dataset__dataroot=dataroot,
        )

    # features are unsupported in the current hosted server
    app_config.update(
        diffexp__enable=False,
        user_annotations__enable=False,
        embeddings__enable_reembedding=False,
        multi_dataset__allowed_matrix_types=["cxg"],
    )

    matrix_data_cache_manager = MatrixDataCacheManager()
    app_config.complete_config(matrix_data_cache_manager, logging.info)
    user_annotations = app_config.user_annotations

    server = Server(matrix_data_cache_manager, user_annotations, app_config)

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
