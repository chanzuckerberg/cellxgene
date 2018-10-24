import os
import warnings

import numpy as np
from pandas import DataFrame, Series
import scanpy.api as sc
from scipy import stats

from server.app.app import cache
from server.app.driver.driver import CXGDriver
from server.app.util.constants import Axis, DEFAULT_TOP_N, DiffExpMode
from server.app.util.utils import FilterError, InteractiveError, PrepareError

"""
Sort order for methods
1. Initialize
2. Helper
3. Filter
4. Data & Metadata
5. Computation
"""


class ScanpyEngine(CXGDriver):

    def __init__(self, data, layout_method=None, diffexp_method=None, obs_names=None, var_names=None):
        super().__init__(data, layout_method=layout_method, diffexp_method=diffexp_method)
        self._alias_annotation_names(obs_names, var_names)
        self._validate_data_types()
        self.cell_count = self.data.shape[0]
        self.gene_count = self.data.shape[1]
        self._create_schema()

    def _alias_annotation_names(self, obs_names, var_names):
        """
        Do all user-specified annotaiton aliasing.
        As a *critical* side-effect, ensure the indices are simple number ranges
        """
        aliased_names = {'obs': obs_names, 'var': var_names}
        for ax in Axis:
            name = aliased_names[ax]
            if name == 'name':
                continue

            ax_name = str(ax)
            df_axis = getattr(self.data, ax_name)
            if name is None:
                # reset index to simple range; alias 'name' to point at the
                # previously specified index.
                df_axis = df_axis.reset_index().rename(columns={'index': 'name'})
            elif name in df_axis.columns:
                if name not in df_axis.columns:
                    raise KeyError(f"Annotation name {name}, specified in --{ax_name}-name does not exist.")
                if not df_axis[name].is_unique:
                    raise KeyError(f"Values in -{ax_name}-name must be unique. Please prepare data to contain unique values.")
                # reset index to simple range; alias user-specified annotation to 'name'
                df_axis = df_axis.reset_index(drop=True).rename(columns={name: 'name'})
            else:
                raise KeyError(f"Annotation name {name}, specified in --{ax_name}_name does not exist.")
            setattr(self.data, ax_name, df_axis)

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
        scanpy_group.add_argument("--obs-names",
                                  help="Annotation name to use as unique, human-readable observation name")
        scanpy_group.add_argument("--var-names", help="Annotation name to use as unique, human-readable variable name")
        scanpy_group.add_argument("data_directory", metavar="dir", help="Directory containing data and schema file")
        scanpy_group.set_defaults(func=invocation_function)
        return scanpy_group

    @staticmethod
    def _load_data(data):
        # See https://scanpy.readthedocs.io/en/latest/api/scanpy.api.read.html
        # Based upon this advice, setting cache=True parameter
        # Note: as of current scanpy/anndata release, setting backed='r' will
        # result in an error.
        return sc.read(os.path.join(data, "data.h5ad"), cache=True)

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

    @staticmethod
    def _nan_to_one(values):
        """
        Replaces NaN values with 1
        :param values: numpy ndarray
        :return: ndarray
        """
        return np.where(np.isnan(values), 1, values)

    def _validate_data_types(self):
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

    @cache.memoize()
    def annotation(self, filter, axis, fields=None):
        """
         Gets annotation value for each observation
        :param filter: filter: dictionary with filter params
        :param axis: string obs or var
        :param fields: list of keys for annotation to return, returns all annotation values if not set.
        :return: dict: names - list of fields in order, data - list of lists or metadata
        [observation ids, val1, val2...]
        """
        try:
            df = self.filter_dataframe(filter)
        except KeyError as e:
            raise FilterError(f"Error parsing filter: {e}") from e
        df_axis = getattr(df, axis)
        if not fields:
            fields = df_axis.columns.tolist()
        annotations = DataFrame(df_axis[fields], index=df_axis.index)
        return {
            "names": fields,
            "data": annotations.reset_index().values.tolist()
        }

    @cache.memoize()
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
        try:
            df = self.filter_dataframe(filter)
        except KeyError as e:
            raise FilterError(f"Error parsing filter: {e}") from e
        var_idx = df.var.index.tolist()
        obs_idx = df.obs.index.tolist()
        values = df.X
        df_shape = df.shape
        if df_shape[0] == 1:
            values = values[None, :]
        elif df_shape[1] == 1:
            values = values[:, None]
        if axis == Axis.OBS:
            expression = DataFrame(values, index=df.obs.index)
            result = {
                "var": var_idx,
                "obs": expression.reset_index().values.tolist()
            }
        else:
            expression = DataFrame(values.T, index=df.var.index)
            result = {
                "obs": obs_idx,
                "var": expression.reset_index().values.tolist(),
            }
        return result

    def diffexp(self, filter1, filter2, top_n=None, interactive_limit=None):
        """
        Computes the top differentially expressed variables between two observation sets. If dataframes
        contain a subset of variables, then statistics for all variables will be returned, otherwise
        only the top N vars will be returned.
        :param filter1: filter: dictionary with filter params for first set of observations
        :param filter2: filter: dictionary with filter params for second set of observations
        :param top_n: Limit results to top N (Top var mode only)
        :param interactive_limit: -- don't compute if total # genes in dataframes are larger than this
        :return: top genes, stats and expression values for variables
        """
        try:
            df1 = self.filter_dataframe(filter1)
        except KeyError as e:
            raise FilterError(f"Error parsing filter for set 1: {e}") from e
        # TODO df2 should be inverse if not filter2 provided
        try:
            df2 = self.filter_dataframe(filter2)
        except KeyError as e:
            raise FilterError(f"Error parsing filter for set 2: {e}") from e
        # If not the same genes, test is wrong!
        if np.any(df1.var.index != df2.var.index):
            raise ValueError("Variables ares not the same in set1 and set2")
        if interactive_limit and df1.shape[0] + df2.shape[0] > interactive_limit:
            raise InteractiveError("Size of set 1 and 2 is too large for interactive computation")
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
        pval = self._nan_to_one(diffexp_result.pvalue)
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

    @cache.memoize()
    def layout(self, filter, interactive_limit=None):
        """
        Computes a n-d layout for cells through dimensionality reduction.
        :param filter: filter: dictionary with filter params
        :param interactive_limit: -- don't compute if total # genes in dataframes are larger than this
        :return:  [cellid, x, y, ...]
        """
        try:
            df = self.filter_dataframe(filter, include_uns=True)
        except KeyError as e:
            raise FilterError(f"Error parsing filter: {e}") from e
        if interactive_limit and len(df.obs.index) > interactive_limit:
            raise InteractiveError("Size data is too large for interactive computation")
        # TODO Filtering cells is fine, but filtering genes does nothing because the neighbors are
        # calculated using the original vars (geneset) and this doesn’t get updated when you use less.
        # Need to recalculate neighbors (long) if user requests new layout filtered by var
        # TODO for MVP we are pushing computation of layout to preprocessing and not allowing re-layout
        # this will probably change after user feedback
        # getattr(sc.tl, self.layout_method)(df, random_state=123)
        try:
            df_layout = df.obsm[f"X_{self.layout_method}"]
        except ValueError as e:
            raise PrepareError(f"Layout has not been calculated using {self.layout_method}, "
                               f"please prepare your datafile and relaunch cellxgene") from e
        normalized_layout = DataFrame((df_layout - df_layout.min()) / (df_layout.max() - df_layout.min()),
                                      index=df.obs.index)
        return {
            "ndims": normalized_layout.shape[1],
            # reset_index gets obs' id into output
            "coordinates": normalized_layout.reset_index().values.tolist()
        }
