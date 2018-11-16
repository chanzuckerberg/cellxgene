from abc import ABCMeta, abstractmethod

"""
Sort order for methods
1. Initialize
2. Helper
3. Filter
4. Data & Metadata
5. Computation
"""


class CXGDriver(metaclass=ABCMeta):

    def __init__(self, data, args):
        self.data = self._load_data(data)
        self.layout_method = args["layout"]
        self.diffexp_method = args["diffexp"]
        self.max_category_items = args["max_category_items"]
        self.diffexp_lfc_cutoff = args["diffexp_lfc_cutoff"]
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
            features["layout"]["obs"] = {"available": True, "interactiveLimit": 50000}
        if self.diffexp_method:
            features["diffexp"] = {"available": True, "interactiveLimit": 50000}
        if self.cluster:
            features["cluster"] = {"available": True, "interactiveLimit": 50000}
        return features

    @staticmethod
    @abstractmethod
    def _load_data(data):
        pass

    @abstractmethod
    def filter_dataframe(self, filter):
        """
        Filter cells from data and return a subset of the data. They can operate on both obs and var dimension with
        indexing and filtering by annotation value. Filters are combined with the and operator.
        See REST specs for info on filter format:
        https://github.com/chanzuckerberg/cellxgene/blob/master/docs/REST_API.md

        :param filter: dictionary with filter params
        :return: View into scanpy object with cells/genes filtered
        """
        pass

    @abstractmethod
    def annotation(self, filter, axis, fields=None):
        """
        Gets annotation value for each observation
        :param filter: filter: dictionary with filter params
        :param axis: string obs or var
        :param fields: list of keys for annotation to return, returns all annotation values if not set.
        :return: dict: names - list of fields in order, data - list of lists or metadata
        [observation ids, val1, val2...]
        """
        pass

    @abstractmethod
    def data_frame(self, filter, axis):
        """
        Retrieves data for each variable for observations in data frame
        :param filter: filter: dictionary with filter params
        :param axis: string obs or var
        :return: {
            "var": list of variable ids,
            "obs": [cellid, var1 expression, var2 expression, ...],
        }
        """
        pass

    @abstractmethod
    def diffexp_topN(self, obsFilter1, obsFilter2, top_n=None, interactive_limit=None):
        """
        Computes the top N differentially expressed variables between two observation sets. If mode
        is "TOP_N", then stats for the top N
        dataframes
        contain a subset of variables, then statistics for all variables will be returned, otherwise
        only the top N vars will be returned.
        :param obsFilter1: filter: dictionary with filter params for first set of observations
        :param obsFilter2: filter: dictionary with filter params for second set of observations
        :param top_n: Limit results to top N (Top var mode only)
        :param interactive_limit: -- don't compute if total # genes in dataframes are larger than this
        :return: top N genes and corresponding stats
        """
        pass

    @abstractmethod
    def layout(self, filter, interactive_limit=None):
        """
        Computes a n-d layout for cells through dimensionality reduction.
        :param filter: filter: dictionary with filter params
        :param interactive_limit: -- don't compute if total # genes in dataframes are larger than this
        :return:  [cellid, x, y, ...]
        """
        pass
