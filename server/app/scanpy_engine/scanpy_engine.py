import warnings

import numpy as np
from pandas.core.dtypes.dtypes import CategoricalDtype
import scanpy as sc
from scipy import sparse

from server.app.driver.driver import CXGDriver
from server.app.util.constants import Axis, DEFAULT_TOP_N
from server.app.util.errors import (
    FilterError,
    JSONEncodingValueError,
    PrepareError,
    ScanpyFileError,
)
from server.app.util.utils import jsonify_scanpy, requires_data
from server.app.scanpy_engine.diffexp import diffexp_ttest
from server.app.util.fbs.matrix import encode_matrix_fbs

"""
Sort order for methods
1. Initialize
2. Helper
3. Filter
4. Data & Metadata
5. Computation
"""


class ScanpyEngine(CXGDriver):
    def __init__(self, data=None, args={}):
        super().__init__(data, args)
        if self.data:
            self._validate_and_initialize()

    def update(self, data=None, args={}):
        super().__init__(data, args)
        if self.data:
            self._validate_and_initialize()

    @staticmethod
    def _get_default_config():
        return {
            "layout": "umap",
            "diffexp": "ttest",
            "max_category_items": 100,
            "obs_names": None,
            "var_names": None,
            "diffexp_lfc_cutoff": 0.01,
        }

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
                raise KeyError(
                    f"Annotation name {name}, specified in --{ax_name}-name does not exist."
                )
            if not df_axis[name].is_unique:
                raise KeyError(
                    f"Values in -{ax_name}-name must be unique. "
                    "Please prepare data to contain unique values."
                )
            # reset index to simple range; alias user-specified annotation to "name"
            df_axis.reset_index(drop=True, inplace=True)
            df_axis.rename(inplace=True, columns={name: "name"})
        else:
            raise KeyError(
                f"Annotation name {name}, specified in --{ax_name}_name does not exist."
            )

    @staticmethod
    def _can_cast_to_float32(ann):
        if ann.dtype.kind == "f":
            if not np.can_cast(ann.dtype, np.float32):
                warnings.warn(
                    f"Annotation {ann.name} will be converted to 32 bit float and may lose precision."
                )
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

    @requires_data
    def _create_schema(self):
        self.schema = {
            "dataframe": {
                "nObs": self.cell_count,
                "nVar": self.gene_count,
                "type": str(self.data.X.dtype),
            },
            "annotations": {"obs": [], "var": []},
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
                    raise TypeError(
                        f"Annotations of type {curr_axis[ann].dtype} are unsupported by cellxgene."
                    )
                self.schema["annotations"][ax].append(ann_schema)

    def _load_data(self, data):
        # Based on benchmarking, cache=True has no impact on perf.
        # Note: as of current scanpy/anndata release, setting backed='r' will
        # result in an error.  https://github.com/theislab/anndata/issues/79
        try:
            self.data = sc.read(data, cache=True)
        except ValueError:
            raise ScanpyFileError(
                "File must be in the .h5ad format. Please read "
                "https://github.com/theislab/scanpy_usage/blob/master/170505_seurat/info_h5ad.md to "
                "learn more about this format. You may be able to convert your file into this format "
                "using `cellxgene prepare`, please run `cellxgene prepare --help` for more "
                "information."
            )
        except MemoryError:
            raise ScanpyFileError("Error while loading file: out of memory, file is too large"
                                  " for memory available")
        except Exception as e:
            raise ScanpyFileError(
                f"Error while loading file: {e}, File must be in the .h5ad format, please check "
                f"that your input and try again."
            )

    @requires_data
    def _validate_and_initialize(self):
        self._alias_annotation_names(Axis.OBS, self.config["obs_names"])
        self._alias_annotation_names(Axis.VAR, self.config["var_names"])
        self._validate_data_types()
        self._validate_data_calculations()
        self.cell_count = self.data.shape[0]
        self.gene_count = self.data.shape[1]
        self._create_schema()

    @requires_data
    def _validate_data_types(self):
        if self.data.X.dtype != "float32":
            warnings.warn(
                f"Scanpy data matrix is in {self.data.X.dtype} format not float32. "
                f"Precision may be truncated."
            )
        for ax in Axis:
            curr_axis = getattr(self.data, str(ax))
            for ann in curr_axis:
                datatype = curr_axis[ann].dtype
                downcast_map = {
                    "int64": "int32",
                    "uint32": "int32",
                    "uint64": "int32",
                    "float64": "float32",
                }
                if datatype in downcast_map:
                    warnings.warn(
                        f"Scanpy annotation {ax}:{ann} is in unsupported format: {datatype}. "
                        f"Data will be downcast to {downcast_map[datatype]}."
                    )
                if isinstance(datatype, CategoricalDtype):
                    category_num = len(curr_axis[ann].dtype.categories)
                    if category_num > 500 and category_num > self.config['max_category_items']:
                        warnings.warn(
                            f"{str(ax).title()} annotation '{ann}' has {category_num} categories, this may be "
                            f"cumbersome or slow to display. We recommend setting the "
                            f"--max-category-items option to 500, this will hide categorical "
                            f"annotations with more than 500 categories in the UI"
                        )

    @requires_data
    def _validate_data_calculations(self):
        layout_key = f"X_{self.config['layout']}"
        try:
            assert layout_key in self.data.obsm_keys()
        except AssertionError:
            raise PrepareError(
                f"Cannot find a field with coordinates for the {self.config['layout']} layout requested. A different"
                f" layout may have been computed. The requested layout must be pre-calculated and saved "
                f"back in the h5ad file. You can run "
                f"`cellxgene prepare --layout {self.config['layout']} <datafile>` "
                f"to solve this problem. "
            )

    @staticmethod
    def _annotation_filter_to_mask(filter, d_axis, count):
        mask = np.ones((count,), dtype=bool)
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
        mask = np.zeros((count,), dtype=bool)
        for i in filter:
            if type(i) == list:
                mask[i[0]: i[1]] = True
            else:
                mask[i] = True
        return mask

    @staticmethod
    def _axis_filter_to_mask(filter, d_axis, count):
        mask = np.ones((count,), dtype=bool)
        if "index" in filter:
            mask = np.logical_and(
                mask, ScanpyEngine._index_filter_to_mask(filter["index"], count)
            )
        if "annotation_value" in filter:
            mask = np.logical_and(
                mask,
                ScanpyEngine._annotation_filter_to_mask(
                    filter["annotation_value"], d_axis, count
                ),
            )
        return mask

    @requires_data
    def _filter_to_mask(self, filter, use_slices=True):
        if use_slices:
            obs_selector = slice(0, self.data.n_obs)
            var_selector = slice(0, self.data.n_vars)
        else:
            obs_selector = None
            var_selector = None

        if filter is not None:
            if Axis.OBS in filter:
                obs_selector = self._axis_filter_to_mask(
                    filter["obs"], self.data.obs, self.data.n_obs
                )
            if Axis.VAR in filter:
                var_selector = self._axis_filter_to_mask(
                    filter["var"], self.data.var, self.data.n_vars
                )
        return obs_selector, var_selector

    @requires_data
    def annotation_to_fbs_matrix(self, axis, fields=None):
        if axis == Axis.OBS:
            df = self.data.obs
        else:
            df = self.data.var
        if fields is not None and len(fields) > 0:
            df = df[fields]
        return encode_matrix_fbs(df, col_idx=df.columns)

    @staticmethod
    def slice_columns(X, var_mask):
        """
        Slice columns from the matrix X, as specified by the mask
        Semantically equivalent to X[:, var_mask], but handles sparse
        matrices in a more performant manner.
        """
        if var_mask is None:    # noop
            return X
        if sparse.issparse(X):  # use tuned getcol/hstack for performance
            indices = np.nonzero(var_mask)[0]
            cols = [X.getcol(i) for i in indices]
            return sparse.hstack(cols, format="csc")
        else:   # else, just use standard slicing, which is fine for dense arrays
            return X[:, var_mask]

    @requires_data
    def data_frame_to_fbs_matrix(self, filter, axis):
        """
        Retrieves data 'X' and returns in a flatbuffer Matrix.
        :param filter: filter: dictionary with filter params
        :param axis: string obs or var
        :return: flatbuffer Matrix

        Caveats:
        * currently only supports access on VAR axis
        * currently only supports filtering on VAR axis
        """
        if axis != Axis.VAR:
            raise ValueError("Only VAR dimension access is supported")
        try:
            obs_selector, var_selector = self._filter_to_mask(filter, use_slices=False)
        except (KeyError, IndexError, TypeError) as e:
            raise FilterError(f"Error parsing filter: {e}") from e
        if obs_selector is not None:
            raise FilterError("filtering on obs unsupported")

        # Currently only handles VAR dimension
        X = self.slice_columns(self.data._X, var_selector)
        return encode_matrix_fbs(X, col_idx=np.nonzero(var_selector)[0], row_idx=None)

    @requires_data
    def diffexp_topN(self, obsFilterA, obsFilterB, top_n=None, interactive_limit=None):
        if Axis.VAR in obsFilterA or Axis.VAR in obsFilterB:
            raise FilterError("Observation filters may not contain vaiable conditions")
        try:
            obs_mask_A = self._axis_filter_to_mask(
                obsFilterA["obs"], self.data.obs, self.data.n_obs
            )
            obs_mask_B = self._axis_filter_to_mask(
                obsFilterB["obs"], self.data.obs, self.data.n_obs
            )
        except (KeyError, IndexError) as e:
            raise FilterError(f"Error parsing filter: {e}") from e
        if top_n is None:
            top_n = DEFAULT_TOP_N
        result = diffexp_ttest(
            self.data, obs_mask_A, obs_mask_B, top_n, self.config['diffexp_lfc_cutoff']
        )
        try:
            return jsonify_scanpy(result)
        except ValueError:
            raise JSONEncodingValueError(
                "Error encoding differential expression to JSON"
            )

    @requires_data
    def layout_to_fbs_matrix(self):
        """
        Return the default 2-D layout for cells as a FBS Matrix.

        Caveats:
        * does not support filtering
        * only returns Matrix in columnar layout
        """
        try:
            full_embedding = self.data.obsm[f"X_{self.config['layout']}"]
            if full_embedding.shape[1] > 2:
                warnings.warn(f"Warning: found {full_embedding.shape[1]} \
                components of embedding. Using the first two for layout display.")
            df_layout = full_embedding[:, :2]
        except ValueError as e:
            raise PrepareError(
                f"Layout has not been calculated using {self.config['layout']}, "
                f"please prepare your datafile and relaunch cellxgene") from e

        normalized_layout = (df_layout - df_layout.min()) / (df_layout.max() - df_layout.min())
        return encode_matrix_fbs(normalized_layout.astype(dtype=np.float32), col_idx=None, row_idx=None)
