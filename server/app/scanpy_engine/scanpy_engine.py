import warnings

import numpy as np
from pandas import DataFrame
from pandas.core.dtypes.dtypes import CategoricalDtype
import scanpy.api as sc
from scipy import sparse

from server.app.driver.driver import CXGDriver
from server.app.util.constants import Axis, DEFAULT_TOP_N
from server.app.util.errors import FilterError, InteractiveError, PrepareError, ScanpyFileError
from server.app.scanpy_engine.diffexp import diffexp_ttest

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

        # TODO: temporary work-arounds
        if args['nan_to_num']:
            self._IEEE754_special_values_workaround()

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

    @staticmethod
    def _can_cast_to_float32(ann):
        if ann.dtype.kind == "f":
            if not np.can_cast(ann.dtype, np.float32):
                warnings.warn(f"Annotation {ann.name} will be converted to 32 bit float and may lose precision.")
            return True
        return False

    @staticmethod
    def _can_cast_to_int32(ann):
        if ann.dtype.kind in ["i", "u"]:
            if np.can_cast(ann.dtype, np.int32):
                return True
            ii32 = np.iinfo(np.int32)
            if ann.min() >= ii32.min and ann.max() <= ii32.max:
                return True
        return False

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
                dtype = curr_axis[ann].dtype
                data_kind = dtype.kind

                if self._can_cast_to_float32(curr_axis[ann]):
                    ann_schema["type"] = "float32"
                elif self._can_cast_to_int32(curr_axis[ann]):
                    ann_schema["type"] = "int32"
                elif dtype == np.bool_:
                    ann_schema["type"] = "boolean"
                elif data_kind == "O" and dtype == "object":
                    ann_schema["type"] = "string"
                elif data_kind == "O" and dtype == "category":
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

    def _IEEE754_special_values_workaround(self):
        """
        TODO: temporary workaround

        Because all floating point data is serialized to JSON, and JSON has no means of representing
        non-finite, floating point special values (NaN, +/-Infinity, etc), we include this temporary
        work-around.

        This will likely be removed in the future, contingent upon improved marshalling.

        Where non-finite floating point is present in obs, var or X:
        * issue a warning to the user that these values will be convert to finite numbers.
        * set NaN to zero, and Infinities to min/max of the element.
        """

        # annotations
        for ax in Axis:
            curr_axis = getattr(self.data, str(ax))
            for ann in curr_axis:
                dtype = curr_axis[ann].dtype
                if dtype.kind == 'f':
                    finite_idx = np.isfinite(curr_axis[ann])
                    if not finite_idx.all():
                        curr_axis.loc[np.isnan(curr_axis[ann]), ann] = 0
                        curr_axis.loc[np.isneginf(curr_axis[ann]), ann] = curr_axis[ann][finite_idx].min()
                        curr_axis.loc[np.isposinf(curr_axis[ann]), ann] = curr_axis[ann][finite_idx].max()
                        warnings.warn(
                            f"{str(ax).title()} annotation '{ann}' contains floating point NaN or Infinities. "
                            f"These will be converted to finite values."
                        )

        # X
        non_finite_X_found = False
        if sparse.issparse(self.data._X):
            coo = self.data._X.tocoo()
            finite_idx = np.isfinite(coo.data)
            if not finite_idx.all():
                non_finite_X_found = True
                coo.data[np.isnan(coo.data)] = 0
                coo.data[np.isneginf(coo.data)] = np.min(coo.data[finite_idx])
                coo.data[np.isposinf(coo.data)] = np.max(coo.data[finite_idx])
                coo.eliminate_zeros()
                _X = coo.asformat(self.data._X.getformat())
                self.data._X = _X
        else:
            _X = self.data._X
            finite_idx = np.isfinite(_X.flat)
            if not finite_idx.all():
                non_finite_X_found = True
                min_X = _X.flat[finite_idx].min()
                max_X = _X.flat[finite_idx].max()
                _X[np.isnan(_X)] = 0
                _X[np.isneginf(_X)] = min_X
                _X[np.isposinf(_X)] = max_X

        if non_finite_X_found:
            warnings.warn(
                "Dataframe X contains floating point NaN or Infinities. "
                "These will be converted to finite values."
            )

    def filter_dataframe(self, filter):
        """
         Filter cells from data and return a subset of the data. They can operate on both obs and var dimension with
         indexing and filtering by annotation value. Filters are combined with the and operator.
         See REST specs for info on filter format:
         # TODO update this link to swagger when it's done
         https://docs.google.com/document/d/1Fxjp1SKtCk7l8QP9-7KAjGXL0eldi_qEnNT0NmlGzXI/edit#heading=h.8qc9q57amldx

        :param filter: dictionary with filter params
        :return: View into scanpy object with cells/genes filtered
        """
        if not filter:
            return self.data
        obs_selector, var_selector = self._filter_to_mask(filter, use_slices=False)
        data = self._slice(self.data, obs_selector, var_selector)
        return data

    @staticmethod
    def _annotation_filter_to_mask(filter, d_axis, count):
        mask = np.ones((count, ), dtype=bool)
        for v in filter:
            if d_axis[v["name"]].dtype.name in ["boolean", "category", "object"]:
                key_idx = np.in1d(getattr(d_axis, v["name"]), v["values"])
                mask = np.logical_and(mask, key_idx)
            else:
                min_ = v.get("min", None)
                max_ = v.get("max", None)
                if min_ is not None:
                    key_idx = (getattr(d_axis, v["name"]) >= min_).ravel()
                    mask = np.logical_and(mask, key_idx)
                if max_ is not None:
                    key_idx = (getattr(d_axis, v["name"]) <= max_).ravel()
                    mask = np.logical_and(mask, key_idx)
        return mask

    @staticmethod
    def _index_filter_to_mask(filter, count):
        mask = np.zeros((count, ), dtype=bool)
        for i in filter:
            if type(i) == list:
                mask[i[0]:i[1]] = True
            else:
                mask[i] = True
        return mask

    @staticmethod
    def _axis_filter_to_mask(filter, d_axis, count):
        mask = np.ones((count, ), dtype=bool)
        if "index" in filter:
            mask = np.logical_and(mask, ScanpyEngine._index_filter_to_mask(filter["index"], count))
        if "annotation_value" in filter:
            mask = np.logical_and(mask,
                                  ScanpyEngine._annotation_filter_to_mask(filter["annotation_value"],
                                                                          d_axis,
                                                                          count))
        return mask

    def _filter_to_mask(self, filter, use_slices=True):
        if use_slices:
            obs_selector = slice(0, self.data.n_obs)
            var_selector = slice(0, self.data.n_vars)
        else:
            obs_selector = None
            var_selector = None

        if filter is not None:
            if Axis.OBS in filter:
                obs_selector = self._axis_filter_to_mask(filter["obs"], self.data.obs, self.data.n_obs)
            if Axis.VAR in filter:
                var_selector = self._axis_filter_to_mask(filter["var"], self.data.var, self.data.n_vars)
        return obs_selector, var_selector

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
            obs_selector, var_selector = self._filter_to_mask(filter)
        except (KeyError, IndexError) as e:
            raise FilterError(f"Error parsing filter: {e}") from e
        if axis == Axis.OBS:
            obs = self.data.obs[obs_selector]
            if not fields:
                fields = obs.columns.tolist()
            result = {
                "names": fields,
                "data": DataFrame(obs[fields]).to_records(index=True).tolist()
            }
        else:
            var = self.data.var[var_selector]
            if not fields:
                fields = var.columns.tolist()
            result = {
                "names": fields,
                "data": DataFrame(var[fields]).to_records(index=True).tolist()
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
            obs_selector, var_selector = self._filter_to_mask(filter)
        except (KeyError, IndexError) as e:
            raise FilterError(f"Error parsing filter: {e}") from e
        _X = self.data._X[obs_selector, var_selector]
        if sparse.issparse(_X):
            _X = _X.toarray()
        var_index_sliced = self.data.var.index[var_selector]
        obs_index_sliced = self.data.obs.index[obs_selector]
        if axis == Axis.OBS:
            result = {
                "var": var_index_sliced.tolist(),
                "obs": DataFrame(_X, index=obs_index_sliced).to_records(index=True).tolist()
            }
        else:
            result = {
                "obs": obs_index_sliced.tolist(),
                "var": DataFrame(_X.T, index=var_index_sliced).to_records(index=True).tolist()
            }
        return result

    def diffexp_topN(self, obsFilterA, obsFilterB, top_n=None, interactive_limit=None):
        if Axis.VAR in obsFilterA or Axis.VAR in obsFilterB:
            raise FilterError("Observation filters may not contain vaiable conditions")
        try:
            obs_mask_A = self._axis_filter_to_mask(obsFilterA["obs"], self.data.obs, self.data.n_obs)
            obs_mask_B = self._axis_filter_to_mask(obsFilterB["obs"], self.data.obs, self.data.n_obs)
        except (KeyError, IndexError) as e:
            raise FilterError(f"Error parsing filter: {e}") from e
        if top_n is None:
            top_n = DEFAULT_TOP_N
        result = diffexp_ttest(self.data, obs_mask_A, obs_mask_B, top_n, self.diffexp_lfc_cutoff)
        return result

    def layout(self, filter, interactive_limit=None):
        """
        Computes a n-d layout for cells through dimensionality reduction.
        :param filter: filter: dictionary with filter params
        :param interactive_limit: -- don't compute if total # genes in dataframes are larger than this
        :return:  [cellid, x, y, ...]
        """
        try:
            df = self.filter_dataframe(filter)
        except (KeyError, IndexError) as e:
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
