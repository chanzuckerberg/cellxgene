from http import HTTPStatus
from flask import make_response, jsonify

from server import __version__ as cellxgene_version
from server.common.data_locator import DataLocator


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

    checks = [
        (config.single_dataset__datapath is not None or config.multi_dataset__dataroot is not None),
        _is_accessible(config.single_dataset__datapath, config),
        _is_accessible(config.multi_dataset__dataroot, config),
    ]
    health["status"] = "pass" if all(checks) else "fail"
    code = HTTPStatus.OK if health["status"] == "pass" else HTTPStatus.BAD_REQUEST
    return make_response(jsonify(health), code, {"Content-Type": "application/health+json"},)
