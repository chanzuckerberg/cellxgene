
class DataAdaptorTypeFactory:
    """Factory class to create a data adaptor type"""

    data_adaptor_types = {}

    @staticmethod
    def register(name, data_adaptor_type):
        DataAdaptorTypeFactory.data_adaptor_types[name] = data_adaptor_type

    @staticmethod
    def get_type(name):
        return DataAdaptorTypeFactory.data_adaptor_types.get(name)

