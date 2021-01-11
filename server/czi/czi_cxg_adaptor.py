from server.data_cxg.cxg_adaptor import CxgAdaptor
from server.data_common.data_adaptor_factory import DataAdaptorTypeFactory


class CZI_CxgAdaptor(CxgAdaptor):
    """A class that combines cxg dataset data with portal database data"""

    # TODO: this is just a stub now.

    # TODO: introduce a function into the base class that maps a path name for a dataset
    # to a location: e.g. input is something like "/e/dataset.cxg" and output is
    # "s3://hosted-cellxgene/actual_dataset.cxg"

    def __init__(self, data_locator, app_config=None, dataset_config=None):
        super().__init__(data_locator, app_config, dataset_config)

    @staticmethod
    def open(data_locator, app_config, dataset_config=None):
        return CZI_CxgAdaptor(data_locator, app_config, dataset_config)

    @staticmethod
    def set_adaptor_params(params):
        """handle adaptor params"""
        # TODO: params may include a URL of the portal backend and other information like that
        pass

    def get_corpora_props(self):
        # TODO:  send a request to the portal backend to retrieve props from the database
        # instead of the cxg files.
        props = super().get_corpora_props()
        return props


DataAdaptorTypeFactory.register("czi_cxg_adaptor", CZI_CxgAdaptor)

