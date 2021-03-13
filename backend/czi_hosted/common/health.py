from http import HTTPStatus
from flask import make_response, jsonify

from backend.czi_hosted import __version__ as cellxgene_version
from backend.czi_hosted.common import DataLocator


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

    checks = False
    server_config = config.server_config
    if config.is_multi_dataset():
        dataroots = [datapath_dict["dataroot"] for datapath_dict in server_config.multi_dataset__dataroot.values()]
        checks = all([_is_accessible(dataroot, server_config) for dataroot in dataroots])
    else:
        checks = _is_accessible(server_config.single_dataset__datapath, server_config)

    health["status"] = "pass" if checks else "fail"
    code = HTTPStatus.OK if health["status"] == "pass" else HTTPStatus.BAD_REQUEST
    response = make_response(jsonify(health), code)
    response.headers["Content-Type"] = "application/health+json"
    return response
