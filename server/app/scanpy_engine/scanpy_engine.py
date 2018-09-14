import os
import warnings

import numpy as np
from pandas import DataFrame, Series
import scanpy.api as sc
from scipy import stats

from server.app.app import cache
from server.app.driver.driver import CXGDriver
from server.app.util.constants import Axis


class ScanpyEngine(CXGDriver):

    def __init__(self, data, layout_method=None, diffexp_method=None):
        super().__init__(data, layout_method=layout_method, diffexp_method=diffexp_method)
        self._validatate_data_types()
        self._add_mandatory_annotations()
        self.cell_count = self.data.shape[0]
        self.gene_count = self.data.shape[1]
        self._create_schema()
        self.layout(self.data)

    def _create_schema(self):
        self.schema = {
            "dataframe": {
                "nObs": self.cell_count,
                "nVar": self.gene_count,
                "type": str(self.data.X.dtype)
            },
            "annotations": {
                "obs": [],
                "var": []
            }
        }
        for ax in Axis:
            curr_axis = getattr(self.data, str(ax))
            for ann in curr_axis:
                ann_schema = {"name": ann}
                data_kind = curr_axis[ann].dtype.kind
                if data_kind == 'f':
                    ann_schema["type"] = "float32"
                elif data_kind in ['i', 'u']:
                    ann_schema["type"] = "int32"
                elif data_kind == "?":
                    ann_schema["type"] = "boolean"
                elif data_kind == "O" and curr_axis[ann].dtype == "object":
                    ann_schema["type"] = "string"
                elif data_kind == "O" and curr_axis[ann].dtype == "category":
                    ann_schema["type"] = "categorical"
                    ann_schema["categories"] = curr_axis[ann].dtype.categories.tolist()
                else:
                    raise TypeError(f"Annotations of type {curr_axis[ann].dtype} are unsupported by cellxgene.")
                self.schema["annotations"][ax].append(ann_schema)

    @classmethod
    def add_to_parser(cls, subparsers, invocation_function):
        scanpy_group = subparsers.add_parser("scanpy", help="run cellxgene using the scanpy engine")
        # TODO these choices should be generated from the actual available methods see GH issue #94
        scanpy_group.add_argument("-l", "--layout", choices=["umap", "tsne"], default="umap",
                                  help="Algorithm to use for graph layout")
        scanpy_group.add_argument("-d", "--diffexp", choices=["ttest"], default="ttest",
                                  help="Algorithm to used to calculate differential expression")
        scanpy_group.add_argument("data_directory", metavar="dir", help="Directory containing data and schema file")
        scanpy_group.set_defaults(func=invocation_function)
        return scanpy_group

    @staticmethod
    def _load_data(data):
        return sc.read(os.path.join(data, "data.h5ad"))

    def _add_mandatory_annotations(self):
        # ensure gene
        self.data.var["name"] = Series(list(self.data.var.index), dtype="unicode_", index=self.data.var.index)
        self.data.var.index = Series(list(range(self.data.var.shape[0])), dtype="category")
        # ensure cell name
        self.data.obs["name"] = Series(list(self.data.obs.index), dtype="unicode_", index=self.data.obs.index)
        self.data.obs.index = Series(list(range(self.data.obs.shape[0])), dtype="category")

    def _validatate_data_types(self):
        if self.data.X.dtype != "float32":
            warnings.warn(f"Scanpy data matrix is in {self.data.X.dtype} format not float32. "
                          f"Precision may be truncated.")
        for ax in Axis:
            curr_axis = getattr(self.data, str(ax))
            for ann in curr_axis:
                datatype = curr_axis[ann].dtype
                downcast_map = {'int64': 'int32',
                                'uint32': 'int32',
                                'uint64': 'int32',
                                'float64': 'float32',
                                }
                if datatype in downcast_map:
                    warnings.warn(f"Scanpy annotation {ax}:{ann} is in unsupported format: {datatype}. "
                                  f"Data will be downcast to {downcast_map[datatype]}.")

    def cells(self):
        return self.data.obs.index.tolist()

    def genes(self):
        return self.data.var.index.tolist()

    # Can't seem to cache a view of a dataframe, need to investigate why
    def filter_dataframe(self, filter):
        """
         Filter cells from data and return a subset of the data. They can operate on both obs and var dimension with
         indexing and filtering by annotation value. Filters are combined with the and operator.
         See REST specs for info on filter format:
         # TODO update this link to swagger when it's done
         https://docs.google.com/document/d/1Fxjp1SKtCk7l8QP9-7KAjGXL0eldi_qEnNT0NmlGzXI/edit#heading=h.8qc9q57amldx

        :param filter: dictionary with filter parames
        :return: View into scanpy object with cells/genes filtered
        """
        cells_idx = np.ones((self.cell_count,), dtype=bool)
        genes_idx = np.ones((self.gene_count,), dtype=bool)
        if Axis.OBS in filter:
            if "index" in filter["obs"]:
                cells_idx = self._filter_index(filter["obs"]["index"], cells_idx, Axis.OBS)
            if "annotation_value" in filter["obs"]:
                cells_idx = self._filter_annotation(filter["obs"]["annotation_value"], cells_idx, Axis.OBS)
        if Axis.VAR in filter:
            if "index" in filter["var"]:
                genes_idx = self._filter_index(filter["var"]["index"], genes_idx, Axis.VAR)
            if "annotation_value" in filter["var"]:
                genes_idx = self._filter_annotation(filter["var"]["annotation_value"], genes_idx, Axis.VAR)
        # Due to anndata issues we can't index into cells and genes at the same time
        cell_data = self.data[cells_idx, :]
        data = cell_data[:, genes_idx]
        data.uns = cell_data.uns
        return data

    def _filter_index(self, filter, index, axis):
        """
        Filter data based on index. ex. [1, 3, [111:200]]
        :param filter: subset of filter dict for obs/var:index
        :param index: np logical vector containing true for passing false for failing filter
        :param axis: Axis
        :return: np logical vector for whether the data passes the filter
        """
        if axis == Axis.OBS:
            count_ = self.cell_count
        elif axis == Axis.VAR:
            count_ = self.gene_count
        idx_filter = np.zeros((count_,), dtype=bool)
        for i in filter:
            if type(i) == list:
                idx_filter[i[0]:i[1]] = True
            else:
                idx_filter[i] = True
        return np.logical_and(index, idx_filter)

    def _filter_annotation(self, filter, index, axis):
        """
        Filter data based on annotation value
        :param filter: subset of filter dict for obs/var:annotation_value
        :param index: np logical vector containing true for passing false for failing filter
        :param axis: string obs or var
        :return: np logical vector for whether the data passes the filter
        """
        d_axis = getattr(self.data, axis.value)
        for v in filter:
            if d_axis[v["name"]].dtype.name in ["boolean", "category", "string"]:
                key_idx = np.in1d(getattr(d_axis, v["name"]), v["values"])
                index = np.logical_and(index, key_idx)
            else:
                min_ = v.get("min", None)
                max_ = v.get("max", None)
                if min_ is not None:
                    key_idx = (getattr(d_axis, v["name"]) >= min_).ravel()
                    index = np.logical_and(index, key_idx)
                if max_ is not None:
                    key_idx = (getattr(d_axis, v["name"]) <= max_).ravel()
                    index = np.logical_and(index, key_idx)
        return index

    @cache.memoize()
    def annotation(self, df, fields=None):
        """
         Gets annotation value for each observation

        :param df: from filter_cells, dataframe
        :param fields: list of keys for annotation to return, returns all annotation values if not set.
        :return: dict: names - list of fields in order, data - list of lists or metadata
        [observation ids, val1, val2...]
        """
        if not fields:
            fields = df.obs.columns.tolist()
        annotations = DataFrame(df.obs[fields], index=df.obs.index)
        return {
            "names": fields,
            "data": annotations.reset_index().values.tolist()
        }

    @cache.memoize()
    def layout(self, df):
        """
        Computes a n-d layout for cells through dimensionality reduction.
        :param df: from filter_cells, dataframe
        :return:  [cellid, x, y, ...]
        """
        # TODO Filtering cells is fine, but filtering genes does nothing because the neighbors are
        # calculated using the original vars (geneset) and this doesn’t get updated when you use less.
        # Need to recalculate neighbors (long) if user requests new layout filtered by var
        getattr(sc.tl, self.layout_method)(df, random_state=123)
        df_layout = df.obsm[f"X_{self.layout_method}"]
        normalized_layout = DataFrame((df_layout - df_layout.min()) / (df_layout.max() - df_layout.min()),
                                      index=df.obs.index)
        return {
            "ndims": normalized_layout.shape[1],
            # reset_index gets obs' id into output
            "coordinates": normalized_layout.reset_index().values.tolist()
        }

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
