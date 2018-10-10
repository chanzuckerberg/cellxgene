import os
import warnings

import numpy as np
from pandas import DataFrame, Series
import scanpy.api as sc
from scipy import stats

# TODO fix memoization so that it correctly identifies the same request
# from server.app.app import cache
from server.app.driver.driver import CXGDriver
from server.app.util.constants import Axis, DEFAULT_TOP_N, DiffExpMode

"""
Sort order for methods
1. Initialize
2. Helper
3. Filter
4. Data & Metadata
5. Computation
"""


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

    @staticmethod
    def _top_sort(values, sort_order, top_n=None):
        """
        Sorts an iterable in sort order limited by top_n
        :param values: iterable of values to sort
        :param sort_order: ndarray order to sort in
        :param top_n: cutoff number to return
        :return: values sorted by sort_order limited by top_n
        """
        return values[sort_order][:top_n]

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
    def filter_dataframe(self, filter, include_uns=False):
        """
         Filter cells from data and return a subset of the data. They can operate on both obs and var dimension with
         indexing and filtering by annotation value. Filters are combined with the and operator.
         See REST specs for info on filter format:
         # TODO update this link to swagger when it's done
         https://docs.google.com/document/d/1Fxjp1SKtCk7l8QP9-7KAjGXL0eldi_qEnNT0NmlGzXI/edit#heading=h.8qc9q57amldx

        :param filter: dictionary with filter params
        :param include_uns: bool, include unstructured annotations
        :return: View into scanpy object with cells/genes filtered
        """
        if not filter:
            return self.data
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
        # TODO: tmp hack to avoid problems with filter that is limited to single gene
        if include_uns:
            data.uns = cell_data.uns
        return data

    def _filter_index(self, filter, index, axis):
        """
        Filter data based on index. ex. [1, 3, [111:200]]
        :param filter: subset of filter dict for obs/var:index
        :param index: np logical vector containing true for passing false for failing filter
        :param axis: string obs or var
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
            if d_axis[v["name"]].dtype.name in ["boolean", "category", "object"]:
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

    # @cache.memoize()
    def annotation(self, df, axis, fields=None):
        """
         Gets annotation value for each observation
        :param df: from filter_cells, dataframe
        :param axis: string obs or var
        :param fields: list of keys for annotation to return, returns all annotation values if not set.
        :return: dict: names - list of fields in order, data - list of lists or metadata
        [observation ids, val1, val2...]
        """
        df_axis = getattr(df, axis)
        if not fields:
            fields = df_axis.columns.tolist()
        annotations = DataFrame(df_axis[fields], index=df_axis.index)
        return {
            "names": fields,
            "data": annotations.reset_index().values.tolist()
        }

    # @cache.memoize()
    def data_frame(self, df, axis):
        """
        Retrieves data for each variable for observations in data frame
        :param df: from filter_cells, dataframe
        :param axis: string obs or var
        :return: {
            "var": list of variable ids,
            "obs": [cellid, var1 expression, var2 expression, ...],
        }
        """
        var_idx = df.var.index.tolist()
        obs_idx = df.obs.index.tolist()
        values = df.X
        df_shape = df.shape
        if df_shape[0] == 1:
            values = values[None, :]
        elif df_shape[1] == 1:
            values = values[:, None]
        if axis == Axis.OBS:
            expression = DataFrame(values, index=obs_idx)
            result = {
                "var": var_idx,
                "obs": expression.reset_index().values.tolist()
            }
        else:
            expression = DataFrame(values.T, index=var_idx)
            result = {
                "obs": obs_idx,
                "var": expression.reset_index().values.tolist(),
            }
        return result

    # @cache.memoize()
    def diffexp(self, df1, df2, top_n=None):
        """
        Computes the top differentially expressed variables between two observation sets. If dataframes
        contain a subset of variables, then statistics for all variables will be returned, otherwise
        only the top N vars will be returned.
        :param df1: from filter_cells, dataframe containing first set of observations
        :param df2: from filter_cells, dataframe containing second set of observations
        :param top_n: Limit results to top N (Top var mode only)
        :return: top genes, stats and expression values for variables
        """
        # If not the same genes, test is wrong!
        if np.any(df1.var.index != df2.var.index):
            raise ValueError("Variables ares not the same in set1 and set2")

        # If not all genes, they used a var filter
        if df1.var.shape[0] < self.gene_count:
            mode = DiffExpMode.VAR_FILTER
            if top_n:
                raise Warning("Top N was specified but will not be used in 'Var Filter' mode")
        else:
            mode = DiffExpMode.TOP_N
            if not top_n:
                top_n = DEFAULT_TOP_N

        genes_idx = df1.var.index
        diffexp_result = stats.ttest_ind(df1.X, df2.X)
        pval = diffexp_result.pvalue
        bonferroni_pval = 1 - (1 - pval) ** self.gene_count
        ave_exp_set1 = np.mean(df1.X, axis=0)
        ave_exp_set2 = np.mean(df2.X, axis=0)
        ave_diff = ave_exp_set1 - ave_exp_set2
        if mode == DiffExpMode.TOP_N:
            sort_order = np.argsort(np.abs(diffexp_result.statistic))[::-1]
            # If top_n > length it will just return length
            genes = self._top_sort(genes_idx, sort_order, top_n)
            pval = self._top_sort(pval, sort_order, top_n)
            bonferroni_pval = self._top_sort(bonferroni_pval, sort_order, top_n)
            ave_exp_set1 = self._top_sort(ave_exp_set1, sort_order, top_n)
            ave_exp_set2 = self._top_sort(ave_exp_set2, sort_order, top_n)
            ave_diff = self._top_sort(ave_diff, sort_order, top_n)

        # varIndex, avgDiff,  pVal, pValAdj, set1AvgExp, set2AvgExp
        result = []
        for i in range(len(genes)):
            result.append([genes[i], ave_diff[i], pval[i], bonferroni_pval[i], ave_exp_set1[i], ave_exp_set2[i]])
        # Results need to be returned in var index order
        return sorted(result, key=lambda gene: gene[0])

    # @cache.memoize()
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
