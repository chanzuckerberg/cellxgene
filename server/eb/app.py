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
    logging.exception("Exception importing server modules")
    sys.exit(1)

try:
    dataroot = os.getenv("CXG_DATAROOT")
    if dataroot is None:
        logging.error("CXG_DATAROOT environment variable must be set")
        sys.exit(1)

    app_config = AppConfig(
        datapath = None,
        dataroot = dataroot,
        title = "",
        about = None,
        scripts = [],
        layout = [],
        max_category_items = 100,
        diffexp_lfc_cutoff = 0.01,
        obs_names = None,
        var_names = None,
        anndata_backed = False,
        disable_diffexp = False)

    matrix_data_cache_manager = MatrixDataCacheManager()
    annotations = None

    server = Server(matrix_data_cache_manager, annotations, app_config)

    debug = False
    application = server.app

except Exception:
    logging.exception("Caught exception during initialization")
    sys.exit(1)

logging.info(f"starting server with CXG_DATAROOT={dataroot}")

if __name__ == "__main__":
    try:
        application.run(debug=debug, threaded=not debug, use_debugger=False)
    except Exception:
        logging.exception("Caught exception during server run")
        sys.exit(1)
