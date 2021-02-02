from enum import Enum
from local_server.common.errors import DatasetAccessError
from local_server.common.data_locator import DataLocator
from http import HTTPStatus


class MatrixDataType(Enum):
    H5AD = "h5ad"
    UNKNOWN = "unknown"


class MatrixDataLoader(object):
    def __init__(self, location, matrix_data_type=None, app_config=None):
        """ location can be a string or DataLocator """
        region_name = None if app_config is None else app_config.server_config.data_locator__s3__region_name
        self.location = DataLocator(location, region_name=region_name)
        if not self.location.exists():
            raise DatasetAccessError("Dataset does not exist.", HTTPStatus.NOT_FOUND)

        # matrix_data_type is an enum value of type MatrixDataType
        self.matrix_data_type = matrix_data_type
        # matrix_type is a DataAdaptor type, which corresonds to the matrix_data_type
        self.matrix_type = None

        if matrix_data_type is None:
            self.matrix_data_type = self.__matrix_data_type()

        if not self.__matrix_data_type_allowed(app_config):
            raise DatasetAccessError("Dataset does not have an allowed type.")

        if self.matrix_data_type == MatrixDataType.H5AD:
            from local_server.data_anndata.anndata_adaptor import AnndataAdaptor

            self.matrix_type = AnndataAdaptor

    def __matrix_data_type(self):
        if self.location.path.endswith(".h5ad"):
            return MatrixDataType.H5AD
        else:
            return MatrixDataType.UNKNOWN

    def __matrix_data_type_allowed(self, app_config):
        return self.matrix_data_type != MatrixDataType.UNKNOWN

    def pre_load_validation(self):
        if self.matrix_data_type == MatrixDataType.UNKNOWN:
            raise DatasetAccessError("Dataset does not have a recognized type: .h5ad")
        self.matrix_type.pre_load_validation(self.location)

    def file_size(self):
        return self.matrix_type.file_size(self.location)

    def open(self, app_config, dataset_config=None):
        # create and return a DataAdaptor object
        return self.matrix_type.open(self.location, app_config, dataset_config)
