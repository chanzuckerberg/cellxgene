from abc import ABCMeta, abstractmethod


class CXGDriver(metaclass=ABCMeta):

    def __init__(self, data, layout_method=None, diffexp_method=None):
        self.data = self._load_data(data)
        self.layout_method = layout_method
        self.diffexp_method = diffexp_method
        self.cluster = None

    @property
    def features(self):
        features = {
            "cluster": {"available": False},
            "layout": {
                "obs": {"available": False},
                "var": {"available": False},
            },
            "diffexp": {"available": False}
        }
        # TODO - Interactive limit should be generated from the actual available methods see GH issue #94
        if self.layout_method:
            # TODO handle "var" when gene layout becomes available
            features["layout"]["obs"] = {"available": True, "interactiveLimit": 15000}
        if self.diffexp_method:
            features["diffexp"] = {"available": True, "interactiveLimit": 5000}
        if self.cluster:
            features["cluster"] = {"available": True, "interactiveLimit": 45000}
        return features

    @staticmethod
    @abstractmethod
    def _load_data(data):
        pass

    @abstractmethod
    def cells(self):
        pass

    @abstractmethod
    def genes(self):
        pass

    @abstractmethod
    def filter_dataframe(self, filter):
        """
         Filter cells from data and return a subset of the data. They can operate on both obs and var dimension with
         indexing and filtering by annotation value. Filters are combined with the and operator.
         See REST specs for info on filter format:
         https://docs.google.com/document/d/1Fxjp1SKtCk7l8QP9-7KAjGXL0eldi_qEnNT0NmlGzXI/edit#heading=h.8qc9q57amldx

        :param filter: dictionary with filter parames
        :return: View into scanpy object with cells/genes filtered
        """
        pass

    @abstractmethod
    def annotation(self, df, axis, fields=None):
        """
         Gets annotation value for each observation

        :param axis:
        :param df: from filter_cells, dataframe
        :param fields: list of keys for annotation to return, returns all annotation values if not set.
        :return: dict: names - list of fields in order, data - list of lists or metadata [idx, val1, val2...]
        """
        pass

    @abstractmethod
    def layout(self, df):
        """
        Computes a n-d layout for cells through dimensionality reduction.
        :param df: from filter_cells, dataframe
        :return:  [cellid, x, y]
        """
        pass

    @abstractmethod
    def diffexp(self, df1, df2):
        """
        Computes the top differentially expressed genes between two clusters
        :param df1: from filter_cells, dataframe containing first set of cells
        :param df2: from filter_cells, dataframe containing second set of cells
        :return: top genes, stats and expression values for top genes
        """
        pass

    @abstractmethod
    def expression(self, df):
        """
        Retrieves expression for each gene for cells in data frame
        :param df: from filter_cells, dataframe
        :return: {
            "var": list of variable ids,
            "obs": [cellid, var1 expression, var2 expression, ...],
        }
        """
        pass
