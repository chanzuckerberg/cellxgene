from abc import ABCMeta, abstractmethod

from server.common.errors import DisabledFeatureError
from server.common.utils.type_conversion_utils import get_schema_type_hint_of_array
from server.common.genesets import write_gene_sets_tidycsv


class Annotations(metaclass=ABCMeta):
    """baseclass for annotations and gene sets"""

    def __init__(self, config={}):
        self.config = config

    def user_annotations_enabled(self):
        return self.config.get("user-annotations", False)

    def gene_sets_save_enabled(self):
        return self.config.get("genesets-save", False)

    def check_user_annotations_enabled(self):
        if not self.user_annotations_enabled():
            raise DisabledFeatureError("User annotations are disabled.")

    def check_gene_sets_save_enabled(self):
        if not self.gene_sets_save_enabled():
            raise DisabledFeatureError("User gene sets save is disabled.")

    def get_schema(self, data_adaptor):
        schema = []
        labels = self.read_labels(data_adaptor)
        if labels is not None and not labels.empty:
            for col in labels.columns:
                col_schema = dict(name=col, writable=True)
                col_schema.update(get_schema_type_hint_of_array(labels[col]))
                schema.append(col_schema)

        return schema

    @abstractmethod
    def set_collection(self, name):
        """set or create a new annotation collection"""
        pass

    @abstractmethod
    def read_labels(self, data_adaptor):
        """Return the labels as a pandas.DataFrame"""
        pass

    @abstractmethod
    def write_labels(self, df, data_adaptor):
        """Write the labels (df) to a persistent storage such that it can later be read"""
        pass

    @abstractmethod
    def read_gene_sets(self, data_adaptor):
        """Return the gene sets from persistent storage"""
        pass

    @abstractmethod
    def write_gene_sets(self, gs, data_adaptor):
        """Write the gene sets (gs) to a persistent storage such that it can later be read"""
        pass

    @abstractmethod
    def update_parameters(self, parameters, data_adaptor):
        """Update configuration parameters that describe information about the annotations feature"""
        pass

    @staticmethod
    def gene_sets_to_csv(genesets):
        """
        Convert the internal gene sets format (returned by read_gene_set) into
        the simple Tidy CSV.
        """
        from io import StringIO

        if isinstance(genesets, dict):
            genesets = genesets.values()

        with StringIO() as sio:
            write_gene_sets_tidycsv(sio, genesets)
            return sio.getvalue()

    @staticmethod
    def gene_sets_to_response(genesets):
        """
        Convert the internal gene sets format (returned by read_gene_set) into
        the dict expected by the JSON REST API
        """
        return list(genesets.values())
