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
    def __init__(self, data=None, args={}):
        self.config = self._get_default_config()
        self.config.update(args)
        if data:
            self._load_data(data)
        else:
            self.data = None

    def update(self, data=None, args={}):
        self.config.update(args)
        if data:
            self._load_data(data)

    @staticmethod
    def _get_default_config():
        return {
            "layout": None,
            "max_category_items": None,
            "diffexp_lfc_cutoff": None,
            "disable_diffexp": False,
            "diffexp_may_be_slow": False
        }

    @property
    def features(self):
        features = {
            "cluster": {"available": False},
            "layout": {"obs": {"available": False}, "var": {"available": False}},
            "diffexp": {"available": True, "interactiveLimit": 50000}
        }
        # TODO - Interactive limit should be generated from the actual available methods see GH issue #94
        if self.config["layout"]:
            # TODO handle "var" when gene layout becomes available
            features["layout"]["obs"] = {"available": True, "interactiveLimit": 50000}
        return features

    @abstractmethod
    def get_schema(self):
        """
        Return current schema
        """
        pass

    @abstractmethod
    def _load_data(self, data_locator):
        pass

    @abstractmethod
    def annotation_to_fbs_matrix(self, axis, field=None):
        """
        Gets annotation value for each observation
        :param axis: string obs or var
        :param fields: list of keys for annotation to return, returns all annotation values if not set.
        :return: flatbuffer: in fbs/matrix.fbs encoding
        """
        pass

    @abstractmethod
    def annotation_put_fbs(self, axis, fbs):
        """
        Put/save FBS as user-defined labels
        """
        pass

    @abstractmethod
    def data_frame_to_fbs_matrix(self, filter, axis):
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
    def layout_to_fbs_matrix(self, filter):
        """ same as layout, except returns a flatbuffer """
        pass
