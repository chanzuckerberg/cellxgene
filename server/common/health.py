from http import HTTPStatus
from flask import make_response, jsonify

from server import __version__ as cellxgene_version
from server.common.data_locator import DataLocator


def _is_accessible(path):
    if path is None:
        return True

    try:
        dl = DataLocator(path)
        return dl.exists()
    except RuntimeError:
        return False


def health_check(config):
    """
    simple health check - return HTTP response.
    See https://tools.ietf.org/id/draft-inadarei-api-health-check-01.html
    """
    health = {"status": None, "version": "1", "releaseID": cellxgene_version}

    checks = [_is_accessible(config.datapath), _is_accessible(config.dataroot)]
    health["status"] = "pass" if all(checks) else "fail"
    code = HTTPStatus.OK if health["status"] == "pass" else HTTPStatus.BAD_REQUEST
    return make_response(jsonify(health), code, {"Content-Type": "application/health+json"},)
