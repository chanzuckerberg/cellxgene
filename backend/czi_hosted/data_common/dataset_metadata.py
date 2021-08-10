import json
import logging
from http import HTTPStatus

import requests
from backend.common.utils.utils import path_join
from backend.common.errors import DatasetAccessError
import backend.czi_hosted.common.rest as common_rest
from backend.czi_hosted.common.config.app_config import AppConfig
from backend.czi_hosted.common.config.server_config import ServerConfig


def request_dataset_metadata_from_data_portal(data_portal_api_base: str, explorer_url: str):
    """
    Check the data portal metadata api for datasets stored under the given url_path
    If present return dataset metadata object else return None
    """
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    try:
        response = requests.get(
            url=f"{data_portal_api_base}/datasets/meta?url={explorer_url}",
            headers=headers
        )
        if response.status_code == 200:
            dataset_identifiers = json.loads(response.body)
            return dataset_identifiers
        else:
            return None
    except Exception:
        return None

def extrapolate_dataset_location_from_config(server_config: ServerConfig, dataset_explorer_location: str):
    """
    Use the dataset_explorer_location and the server config to determine where the dataset is stored
    """
    # TODO @mdunitz remove after fork, update config to remove single_dataset option, the multiroot lookup will need to
    #  remain while we support covid 19 cell atlas
    if server_config.single_dataset__datapath:
        datapath = server_config.single_dataset__datapath
        return datapath
    else:
        dataset_explorer_location = dataset_explorer_location.split("/")
        dataset = dataset_explorer_location.pop(-1)
        url_dataroot = "/".join(dataset_explorer_location)
        dataroot = None
        for key, dataroot_dict in server_config.multi_dataset__dataroot.items():
            if dataroot_dict["base_url"] == url_dataroot:
                dataroot = dataroot_dict["dataroot"]
                break
        if dataroot is None:
            raise DatasetAccessError(f"Invalid dataset {url_dataroot}/{dataset}")
        datapath = path_join(dataroot, dataset)
        # path_join returns a normalized path.  Therefore it is
        # sufficient to check that the datapath starts with the
        # dataroot to determine that the datapath is under the dataroot.
        if not datapath.startswith(dataroot):
            raise DatasetAccessError(f"Invalid dataset {url_dataroot}/{dataset}")
        return datapath


def get_dataset_metadata_for_explorer_location(dataset_explorer_location: str, app_config: AppConfig):
    """
    Given the dataset access location(the path to view the dataset in the explorer including the dataset root,
    also used as the cache key) and the explorer web base url (from the app_config) this function returns the location
    of the dataset, along with additional metadata.
    The dataset location is is either retrieved from the data portal, or built based on the dataroot information stored
    in the server config.
    In the case of a single dataset the dataset location is pulled directly from the server_config.
    """
    if app_config.server_config.data_locator__api_base:
        explorer_url_path = f"{app_config.server_config.get_web_base_url()}/{dataset_explorer_location}"
        dataset_metadata = request_dataset_metadata_from_data_portal(
            data_portal_api_base=app_config.server_config.data_locator__api_base,
            explorer_url=explorer_url_path
        )
        if dataset_metadata:
            return dataset_metadata
    server_config = app_config.server_config
    dataset_metadata = {
        "collection_id": None,
        "collection_visibility": None,
        "dataset_id": None,
        "s3_uri": None,
        "tombstoned": False
    }
    dataset_metadata["s3_uri"] = extrapolate_dataset_location_from_config(server_config=server_config, dataset_explorer_location=dataset_explorer_location)
    if dataset_metadata["s3_uri"] is None:
        return common_rest.abort_and_log(HTTPStatus.BAD_REQUEST, "Invalid dataset NONE", loglevel=logging.INFO)

    return dataset_metadata

