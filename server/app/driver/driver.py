from abc import ABCMeta, abstractmethod


class CXGDriver(metaclass=ABCMeta):
    def __init__(self, data, graph_method=None, diffexp_method=None):
        self.data = self._load_data(data)
        self._set_features(graph_method=graph_method, diffexp_method=diffexp_method, cluster=None)

    def _set_features(self, graph_method=None, diffexp_method=None, cluster=None):
        self.features = {
            "cluster": {"available": False},
            "layout": {
                "obs": {"available": False},
                "var": {"available": False},
            },
            "diffexp": {"available": False}
        }
        # TODO - Interactive limit should be generated from the actual available methods see GH issue #94
        if graph_method:
            # TODO handle "var" when gene layout becomes available
            self.features["layout"]["obs"] = {"available": True, "interactiveLimit": 15000}
        if diffexp_method:
            self.features["diffexp"] = {"available": True, "interactiveLimit": 5000}
        if cluster:
            self.features["cluster"] = {"available": True, "interactiveLimit": 45000}

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
    def metadata(self, df, fields=None):
        """
         Gets metadata key:value for each cells

        :param df: from filter_cells, dataframe
        :param fields: list of keys for metadata to return, returns all metadata values if not set.
        :return: list of metadata values
        """
        pass

    @abstractmethod
    def create_graph(self, df):
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
        :param df:
        :return: {
            "genes": list of genes,
            "cells": list of cells and expression list,
            "nonzero_gene_count": number of nonzero genes
        }
        """
        pass
