import os

import numpy as np
import scanpy.api as sc
from scipy import stats

from server.app.app import cache
from server.app.driver.driver import CXGDriver
from server.app.util.schema_parse import parse_schema


class ScanpyEngine(CXGDriver):

    def __init__(self, data, schema=None, graph_method="umap", diffexp_method="ttest"):
        self.data = self._load_data(data)
        self.schema = self._load_or_infer_schema(data, schema)
        self._set_cell_names()
        self.cell_count = self.data.shape[0]
        self.gene_count = self.data.shape[1]
        self.graph_method = graph_method
        self.diffexp_method = diffexp_method

    def _set_cell_names(self):
        self.data.obs["cell_name"] = list(self.data.obs.index)

    @classmethod
    def add_to_parser(cls, subparsers, invocation_function):
        scanpy_group = subparsers.add_parser("scanpy", help="run cellxgene using the scanpy engine")
        # TODO these choices should be generated from the actual available methods
        scanpy_group.add_argument("-l", "--layout", choices=["umap", "tsne"], default="umap",
                                  help="Algorithm to use for graph layout")
        scanpy_group.add_argument("-d", "--diffexp", choices=["ttest"], default="ttest",
                                  help="Algorithm to use to calculate differential expression")
        scanpy_group.add_argument("data_directory", metavar="dir", help="Directory containing data and schema file")
        scanpy_group.set_defaults(func=invocation_function)
        return scanpy_group

    @staticmethod
    def _load_data(data):
        return sc.read(os.path.join(data, "data.h5ad"))

    def _load_or_infer_schema(self, data, schema):
        if not os.path.isfile(os.path.join(data, schema)):
            # Initialize with cell name which is built off the index
            data_schema = {
                "CellName": {
                    "type": "string",
                    "variabletype": "categorical",
                    "displayname": "Name",
                    "include": True
                }
            }
            metadata_fields = list(self.data.obs)
            for m in metadata_fields:
                # Since there are many type of float/int in numpy datatypes the kind attribute of a datatype object
                # offers a decent insight into whether it can be lumped in with floats or ints, which is what we
                # care about here.
                data_kind = self.data.obs[m].dtype.kind
                variable_type = "categorical"
                data_type = "string"
                if data_kind == 'f':
                    variable_type = "continuous"
                    data_type = "float"
                elif data_kind in ['i', 'u']:
                    data_type = "int"
                    if self.data.obs[m].nunique() > 50:
                        variable_type = "continuous"
                data_schema[m] = {
                    "type": data_type,
                    "variabletype": variable_type,
                    "displayname": m,
                    "include": True
                }
        else:
            data_schema = parse_schema(os.path.join(data, schema))
        return data_schema

    def cells(self):
        return list(self.data.obs.index)

    def genes(self):
        return self.data.var.index.tolist()

    # Can't seem to cache a view of a dataframe, need to investigate why
    def filter_cells(self, filter):
        """
        Filter cells from data and return a subset of the data
        A filter is a dictionary where the key is a metadatata category
        Value is dictionary
            value_type: int, float, string
            variable_type: continuous, categorical
            query: filter value, for categorical [val1, val2], for continuous {min: x, max:y}
        Filters are combined with the and operator
        :param filter:
        :return: filtered dataframe
        """
        cell_idx = np.ones((self.cell_count,), dtype=bool)
        for key, value in filter.items():
            if value["variable_type"] == "categorical":
                key_idx = np.in1d(getattr(self.data.obs, key), value["query"])
                cell_idx = np.logical_and(cell_idx, key_idx)
            else:
                min_ = value["query"]["min"]
                max_ = value["query"]["max"]
                if min_ is not None:
                    key_idx = np.array((getattr(self.data.obs, key) >= min_).data)
                    cell_idx = np.logical_and(cell_idx, key_idx)
                if max_ is not None:
                    key_idx = np.array((getattr(self.data.obs, key) <= max_).data)
                    cell_idx = np.logical_and(cell_idx, key_idx)
        return self.data[cell_idx, :]

    @cache.memoize()
    def metadata_ranges(self, df=None):
        metadata_ranges = {}
        if not df:
            df = self.data
        for field in self.schema:
            if self.schema[field]["variabletype"] == "categorical":
                group_by = field
                if group_by == "CellName":
                    group_by = "cell_name"
                metadata_ranges[field] = {"options": df.obs.groupby(group_by).size().to_dict()}
            else:
                metadata_ranges[field] = {
                    "range": {
                        "min": df.obs[field].min(),
                        "max": df.obs[field].max()
                    }
                }
        return metadata_ranges

    @cache.memoize()
    def metadata(self, df, fields=None):
        """
         Gets metadata key:value for each cells

        :param df: from filter_cells, dataframe
        :param fields: list of keys for metadata to return, returns all metadata values if not set.
        :return: list of metadata values
        """
        metadata = df.obs.to_dict(orient="records")
        for idx in range(len(metadata)):
            metadata[idx]["CellName"] = metadata[idx].pop("cell_name", None)
        return metadata

    @cache.memoize()
    def create_graph(self, df):
        """
        Computes a n-d layout for cells through dimensionality reduction.
        :param df: from filter_cells, dataframe
        :return:  [cellid, x, y]
        """
        getattr(sc.tl, self.graph_method)(df, random_state=123)
        graph = df.obsm["X_{graph_method}".format(graph_method=self.graph_method)]
        normalized_graph = (graph - graph.min()) / (graph.max() - graph.min())
        return np.hstack((df.obs["cell_name"].values.reshape(len(df.obs.index), 1), normalized_graph)).tolist()

    @cache.memoize()
    def diffexp(self, cell_list_1, cell_list_2, pval, num_genes):
        """
        Computes the top differentially expressed genes between two clusters
        :param df1: from filter_cells, dataframe containing first set of cells
        :param df2: from filter_cells, dataframe containing second set of cells
        :return: top genes, stats and expression values for top genes
        """
        cells_idx_1 = np.in1d(self.data.obs["cell_name"], cell_list_1)
        cells_idx_2 = np.in1d(self.data.obs["cell_name"], cell_list_2)
        expression_1 = self.data.X[cells_idx_1, :]
        expression_2 = self.data.X[cells_idx_2, :]
        diff_exp = stats.ttest_ind(expression_1, expression_2)
        # TODO break this up into functions
        set1 = np.logical_and(diff_exp.pvalue < pval, diff_exp.statistic > 0)
        set2 = np.logical_and(diff_exp.pvalue < pval, diff_exp.statistic < 0)
        stat1 = diff_exp.statistic[set1]
        stat2 = diff_exp.statistic[set2]
        sort_set1 = np.argsort(stat1)[::-1]
        sort_set2 = np.argsort(stat2)
        pval1 = diff_exp.pvalue[set1][sort_set1]
        pval2 = diff_exp.pvalue[set2][sort_set2]
        mean_ex1_set1 = np.mean(expression_1[:, set1], axis=0)[sort_set1]
        mean_ex2_set1 = np.mean(expression_2[:, set1], axis=0)[sort_set1]
        mean_ex1_set2 = np.mean(expression_1[:, set2], axis=0)[sort_set2]
        mean_ex2_set2 = np.mean(expression_2[:, set2], axis=0)[sort_set2]
        mean_diff1 = mean_ex1_set1 - mean_ex2_set1
        mean_diff2 = mean_ex1_set2 - mean_ex2_set2
        genes_cellset_1 = self.data.var_names[set1][sort_set1]
        genes_cellset_2 = self.data.var_names[set2][sort_set2]
        return {
            "celllist1": {
                "topgenes": genes_cellset_1.tolist()[:num_genes],
                "mean_expression_cellset1": mean_ex1_set1.tolist()[:num_genes],
                "mean_expression_cellset2": mean_ex2_set1.tolist()[:num_genes],
                "pval": pval1.tolist()[:num_genes],
                "ave_diff": mean_diff1.tolist()[:num_genes]
            },
            "celllist2": {
                "topgenes": genes_cellset_2.tolist()[:num_genes],
                "mean_expression_cellset1": mean_ex1_set2.tolist()[:num_genes],
                "mean_expression_cellset2": mean_ex2_set2.tolist()[:num_genes],
                "pval": pval2.tolist()[:num_genes],
                "ave_diff": mean_diff2.tolist()[:num_genes]
            },
        }

    @cache.memoize()
    def expression(self, cells=None, genes=None):
        """
        Retrieves expression for each gene for cells in data frame
        :param df:
        :return: {
            "genes": list of genes,
            "cells": list of cells and expression list,
            "nonzero_gene_count": number of nonzero genes
        }
        """
        if cells:
            cells_idx = np.in1d(self.data.obs["cell_name"], cells)
        else:
            cells_idx = np.ones((self.cell_count,), dtype=bool)
        if genes:
            genes_idx = np.in1d(self.data.var_names, genes)
        else:
            genes_idx = np.ones((self.gene_count,), dtype=bool)
        index = np.ix_(cells_idx, genes_idx)
        expression = self.data.X[index]

        if not genes:
            genes = self.data.var.index.tolist()
        if not cells:
            cells = self.data.obs["cell_name"].tolist()

        cell_data = []
        for idx, cell in enumerate(cells):
            cell_data.append({
                "cellname": cell,
                "e": list(expression[idx]),
            })

        return {
            "genes": genes,
            "cells": cell_data,
            "nonzero_gene_count": int(np.sum(expression.any(axis=0)))
        }
