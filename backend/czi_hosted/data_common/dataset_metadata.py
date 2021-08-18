import json
import logging
from http import HTTPStatus

import requests
from backend.common.utils.utils import path_join
from backend.common.errors import DatasetAccessError
import backend.czi_hosted.common.rest as common_rest


def get_dataset_metadata_from_data_portal(explorer_url: str):
    """
    Check the data portal metadata api for datasets stored under the given url_path
    If present return dataset metadata object else return None
    """
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    try:
        response = requests.get(
            url=explorer_url,
            headers=headers
        )
        if response.status_code == 200:
            dataset_identifiers = json.loads(response.body)
            return dataset_identifiers
        else:
            return None
    except Exception:
        return None


def get_dataset_metadata(location: str, **kwargs):
    """
    Given the dataset access location(the path to view the dataset in the explorer including the dataset root,
    also used as the cache key) and the explorer web base url (in the app_config)
     return a dataset_metadata object with the dataset storage location available under s3_uri
    """
    app_config = kwargs["app_config"]
    explorer_url_path = f"{app_config.server_config.get_web_base_url()}/{location}"
    dataset_metadata = get_dataset_metadata_from_data_portal(explorer_url=explorer_url_path)
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
    # TODO @mdunitz remove after fork, update config to remove single_dataset option, the multiroot lookup will need to remain while we support covid 19 cell atlas
    if location is None:
        datapath = server_config.single_dataset__datapath
        dataset_metadata["s3_uri"] = datapath
    else:
        url_dataroot = location.split("/")[0] # TODO check that this returns dataroot (called on e/dataset_id.cxg not /e/dataset_id.cxg)
        dataset = location.split("/")[1]
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

        dataset_metadata["s3_uri"] = datapath
    if dataset_metadata["s3_uri"] is None:
        return common_rest.abort_and_log(HTTPStatus.BAD_REQUEST, "Invalid dataset NONE", loglevel=logging.INFO)

    print(dataset_metadata)
    return dataset_metadata

