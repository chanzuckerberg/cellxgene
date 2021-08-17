from enum import Enum

from backend.common.utils.data_locator import DataLocator
from backend.common.errors import DatasetAccessError
from http import HTTPStatus


class MatrixDataType(Enum):
    H5AD = "h5ad"
    CXG = "cxg"
    UNKNOWN = "unknown"


class MatrixDataLoader(object):
    def __init__(self, location, app_config=None, url_dataroot=None, matrix_data_type=None):
        """ location can be a string or DataLocator """
        self.app_config = app_config
        self.url_dataroot = url_dataroot
        self.dataset_config = self.app_config.default_dataset_config
        if self.url_dataroot:
            dataroot = None
            for key, dataroot_dict in self.app_config.server_config.multi_dataset__dataroot.items():
                if dataroot_dict["base_url"] == self.url_dataroot:
                    dataroot = key
                    break
            if dataroot:
                self.dataset_config = self.app_config.get_dataset_config(dataroot)
        if self.dataset_config is None:
            raise DatasetAccessError("Missing dataset config", HTTPStatus.NOT_FOUND)
        region_name = None if app_config is None else app_config.server_config.data_locator__s3__region_name
        self.location = DataLocator(location, region_name=region_name)
        if not self.location.exists():
            raise DatasetAccessError("Dataset does not exist.", HTTPStatus.NOT_FOUND)

        self.matrix_data_type = matrix_data_type
        # matrix_type is a DataAdaptor type, which corresponds to the matrix_data_type
        self.matrix_type = None
        if self.matrix_data_type is None:
            self.matrix_data_type = self.__matrix_data_type()

        if not self.__matrix_data_type_allowed(app_config):
            raise DatasetAccessError("Dataset does not have an allowed type.")

        if self.matrix_data_type == MatrixDataType.H5AD:
            from backend.czi_hosted.data_anndata.anndata_adaptor import AnndataAdaptor

            self.matrix_type = AnndataAdaptor
        elif self.matrix_data_type == MatrixDataType.CXG:
            from backend.czi_hosted.data_cxg.cxg_adaptor import CxgAdaptor

            self.matrix_type = CxgAdaptor

    # TODO @mdunitz remove when removing conversion code, also remove server_config.multi_dataset__allowed_matrix_types
    def __matrix_data_type(self):
        if self.location.path.endswith(".h5ad"):
            return MatrixDataType.H5AD
        elif ".cxg" in self.location.path:
            return MatrixDataType.CXG
        else:
            return MatrixDataType.UNKNOWN

    def __matrix_data_type_allowed(self, app_config):
        if self.matrix_data_type == MatrixDataType.UNKNOWN:
            return False

        if not app_config:
            return True
        if not app_config.is_multi_dataset():
            return True
        if len(app_config.server_config.multi_dataset__allowed_matrix_types) == 0:
            return True

        for val in app_config.server_config.multi_dataset__allowed_matrix_types:
            try:
                if self.matrix_data_type == MatrixDataType(val):
                    return True
            except ValueError:
                # Check case where multi_dataset_allowed_matrix_type does not have a
                # valid MatrixDataType value.  TODO:  Add a feature to check
                # the AppConfig for errors on startup
                return False

        return False

    def pre_load_validation(self):
        if self.matrix_data_type == MatrixDataType.UNKNOWN:
            raise DatasetAccessError("Dataset does not have a recognized type: .h5ad or .cxg")
        self.matrix_type.pre_load_validation(self.location)

    def file_size(self):
        return self.matrix_type.file_size(self.location)

    def open(self):
        # create and return a DataAdaptor object

        return self.matrix_type.open(self.location, self.app_config, self.dataset_config)

    def validate_and_open(self, dataset_location=None, **kwargs):
        # create and return a DataAdaptor object
        self.pre_load_validation()
        return self.open()
