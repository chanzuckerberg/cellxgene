from http import HTTPStatus
from flask import make_response, jsonify

from backend.server import __version__ as cellxgene_version
from backend.common_utils.data_locator import DataLocator


def _is_accessible(path, config):
    if path is None:
        return True

    try:
        dl = DataLocator(path, region_name=config.data_locator__s3__region_name)
        return dl.exists()
    except RuntimeError:
        return False


def health_check(config):
    """
    simple health check - return HTTP response.
    See https://tools.ietf.org/id/draft-inadarei-api-health-check-01.html
    """
    health = {"status": None, "version": "1", "releaseID": cellxgene_version}

    server_config = config.server_config
    check = _is_accessible(server_config.single_dataset__datapath, server_config)

    health["status"] = "pass" if check else "fail"
    code = HTTPStatus.OK if health["status"] == "pass" else HTTPStatus.BAD_REQUEST
    response = make_response(jsonify(health), code)
    response.headers["Content-Type"] = "application/health+json"
    return response
