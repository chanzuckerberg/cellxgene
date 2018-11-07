import warnings

import numpy as np
from pandas import DataFrame
from pandas.core.dtypes.dtypes import CategoricalDtype
import scanpy.api as sc
from scipy import stats, sparse

from server.app.driver.driver import CXGDriver
from server.app.util.constants import Axis, DEFAULT_TOP_N, DiffExpMode
from server.app.util.errors import FilterError, InteractiveError, PrepareError, ScanpyFileError

"""
Sort order for methods
1. Initialize
2. Helper
3. Filter
4. Data & Metadata
5. Computation
"""


class ScanpyEngine(CXGDriver):

    def __init__(self, data, args):
        super().__init__(data, args)
        self._alias_annotation_names(Axis.OBS, args["obs_names"])
        self._alias_annotation_names(Axis.VAR, args["var_names"])
        self._validate_data_types()
        self._validate_data_calculations()
        self.cell_count = self.data.shape[0]
        self.gene_count = self.data.shape[1]
        self.layout_options = ["umap", "tsne"]
        self.diffexp_options = ["ttest"]
        self._create_schema()

    def _alias_annotation_names(self, axis, name):
        """
        Do all user-specified annotation aliasing.

        As a *critical* side-effect, ensure the indices are simple number ranges
        (accomplished by calling pandas.DataFrame.reset_index())
        """
        if name == "name":
            # a noop, so skip it
            return

        ax_name = str(axis)
        df_axis = getattr(self.data, ax_name)
        if name is None:
            # reset index to simple range; alias "name" to point at the
            # previously specified index.
            df_axis.reset_index(inplace=True)
            df_axis.rename(inplace=True, columns={"index": "name"})
        elif name in df_axis.columns:
            if name not in df_axis.columns:
                raise KeyError(f"Annotation name {name}, specified in --{ax_name}-name does not exist.")
            if not df_axis[name].is_unique:
                raise KeyError(f"Values in -{ax_name}-name must be unique. "
                               "Please prepare data to contain unique values.")
            # reset index to simple range; alias user-specified annotation to "name"
            df_axis.reset_index(drop=True, inplace=True)
            df_axis.rename(inplace=True, columns={name: "name"})
        else:
            raise KeyError(f"Annotation name {name}, specified in --{ax_name}_name does not exist.")

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
                if data_kind == "f":
                    ann_schema["type"] = "float32"
                elif data_kind in ["i", "u"]:
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

    @staticmethod
    def _load_data(data):
        # Based on benchmarking, cache=True has no impact on perf.
        # Note: as of current scanpy/anndata release, setting backed='r' will
        # result in an error.  https://github.com/theislab/anndata/issues/79
        try:
            result = sc.read(data, cache=True)
        except ValueError:
            raise ScanpyFileError("File must be in the .h5ad format. Please read "
                                  "https://github.com/theislab/scanpy_usage/blob/master/170505_seurat/info_h5ad.md to "
                                  "learn more about this format. You may be able to convert your file into this format "
                                  "using `cellxgene prepare`, please run `cellxgene prepare --help` for more "
                                  "information.")
        except Exception as e:
            raise ScanpyFileError(f"Error while loading file: {e}, File must be in the .h5ad format, please check "
                                  f"that your input and try again.")
        return result

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

    @staticmethod
    def _nan_to_zero(values):
        """
        Replaces NaN values with 0
        :param values: numpy ndarray
        :return: ndarray
        """
        return np.where(np.isnan(values), 0, values)

    def _validate_data_types(self):
        if self.data.X.dtype != "float32":
            warnings.warn(f"Scanpy data matrix is in {self.data.X.dtype} format not float32. "
                          f"Precision may be truncated.")
        for ax in Axis:
            curr_axis = getattr(self.data, str(ax))
            for ann in curr_axis:
                datatype = curr_axis[ann].dtype
                downcast_map = {"int64": "int32",
                                "uint32": "int32",
                                "uint64": "int32",
                                "float64": "float32",
                                }
                if datatype in downcast_map:
                    warnings.warn(f"Scanpy annotation {ax}:{ann} is in unsupported format: {datatype}. "
                                  f"Data will be downcast to {downcast_map[datatype]}.")
                if isinstance(datatype, CategoricalDtype):
                    category_num = len(curr_axis[ann].dtype.categories)
                    if category_num > 500 and category_num > self.max_category_items:
                        warnings.warn(
                            f"{str(ax).title()} annotation '{ann}' has {category_num} categories, this may be "
                            f"cumbersome or slow to display. We recommend setting the "
                            f"--max-category-items option to 500, this will hide categorical "
                            f"annotations with more than 500 categories in the UI")

    def _validate_data_calculations(self):
        layout_key = f"X_{self.layout_method}"
        try:
            assert layout_key in self.data.obsm_keys()
        except AssertionError:
            raise PrepareError(
                f"Cannot find a field with coordinates for the {self.layout_method} layout requested. A different"
                f" layout may have been computed. The requested layout must be pre-calculated and saved "
                f"back in the h5ad file. You can run "
                f"`cellxgene prepare --layout {self.layout_method} <datafile>` "
                f"to solve this problem. ")

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

        data = self._slice(self.data, cells_idx, genes_idx)
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

    @staticmethod
    def _slice(data, obs_selector=None, vars_selector=None):
        """
        Slice date using any selector that the AnnData object
        supprots for slicing.  If selector is None, will not slice
        on that axis.

        This method exists to optimize filtering/slicing sparse data that has
        access patterns which impact slicing performance.

        https://docs.scipy.org/doc/scipy/reference/sparse.html
        """
        prefer_row_access = sparse.isspmatrix_csr(data._X) or sparse.isspmatrix_lil(data._X) \
            or sparse.isspmatrix_bsr(data._X)
        if prefer_row_access:
            # Row-major slicing
            if obs_selector is not None:
                data = data[obs_selector, :]
            if vars_selector is not None:
                data = data[:, vars_selector]
        else:
            # Col-major slicing
            if vars_selector is not None:
                data = data[:, vars_selector]
            if obs_selector is not None:
                data = data[obs_selector, :]

        return data

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
        result = {
            "names": fields,
            "data": DataFrame(df_axis[fields]).to_records(index=True).tolist()
        }
        return result

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
            slice = self.filter_dataframe(filter)
        except KeyError as e:
            raise FilterError(f"Error parsing filter: {e}") from e
        # convert sparse slice to dense
        X = slice._X.toarray() if sparse.issparse(slice._X) else slice._X
        if axis == Axis.OBS:
            result = {
                "var": slice.var.index.tolist(),
                "obs": DataFrame(X, index=slice.obs.index).to_records(index=True).tolist()
            }
        else:
            result = {
                "obs": slice.obs.index.tolist(),
                "var": DataFrame(X.T, index=slice.var.index).to_records(index=True).tolist()
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
        # ensure we are using a dense ndarray
        X1 = df1._X.toarray() if sparse.issparse(df1._X) else df1._X
        X2 = df2._X.toarray() if sparse.issparse(df2._X) else df2._X
        diffexp_result = stats.ttest_ind(X1, X2)
        tstats = self._nan_to_zero(diffexp_result.statistic)
        pval = self._nan_to_one(diffexp_result.pvalue)
        bonferroni_pval = 1 - (1 - pval) ** self.gene_count
        ave_exp_set1 = np.mean(X1, axis=0)
        ave_exp_set2 = np.mean(X2, axis=0)
        ave_diff = ave_exp_set1 - ave_exp_set2
        if mode == DiffExpMode.TOP_N:
            sort_order = np.argsort(np.abs(tstats))[::-1]
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
            "coordinates": normalized_layout.to_records(index=True).tolist()
        }
