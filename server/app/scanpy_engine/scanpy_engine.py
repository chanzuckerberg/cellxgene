import warnings

import numpy as np
import pandas
from pandas.core.dtypes.dtypes import CategoricalDtype
import anndata
from scipy import sparse

from server.app.driver.driver import CXGDriver
from server.app.util.constants import Axis, DEFAULT_TOP_N, MAX_LAYOUTS
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
            "layout": [],
            "max_category_items": 100,
            "obs_names": None,
            "var_names": None,
            "diffexp_lfc_cutoff": 0.01,
        }

    @staticmethod
    def _create_unique_column_name(df, col_name_prefix):
        """ given the columns of a dataframe, and a name prefix, return a column name which
            does not exist in the dataframe, AND which is prefixed by `prefix`

            The approach is to append a numeric suffix, starting at zero and increasing by
            one, until an unused name is found (eg, prefix_0, prefix_1, ...).
        """
        suffix = 0
        while f"{col_name_prefix}{suffix}" in df:
            suffix += 1
        return f"{col_name_prefix}{suffix}"

    def _alias_annotation_names(self):
        """
        The front-end relies on the existance of a unique, human-readable
        index for obs & var (eg, var is typically gene name, obs the cell name).
        The user can specify these via the --obs-names and --var-names config.
        If they are not specified, use the existing index to create them, giving
        the resulting column a unique name (eg, "name").

        In both cases, enforce that the result is unique, and communicate the
        index column name to the front-end via the obs_names and var_names config
        (which is incorporated into the schema).
        """
        for (ax_name, config_name) in ((Axis.OBS, "obs_names"), (Axis.VAR, "var_names")):
            name = self.config[config_name]
            df_axis = getattr(self.data, str(ax_name))
            if name is None:
                # Default: create unique names from index
                if not df_axis.index.is_unique:
                    raise KeyError(
                        f"Values in {ax_name}.index must be unique. "
                        "Please prepare data to contain unique index values, or specify an "
                        "alternative with --{ax_name}-name."
                    )
                name = self._create_unique_column_name(df_axis.columns, "name_")
                self.config[config_name] = name
                # reset index to simple range; alias name to point at the
                # previously specified index.
                df_axis.rename_axis(name, inplace=True)
                df_axis.reset_index(inplace=True)
            elif name in df_axis.columns:
                # User has specified alternative column for unique names, and it exists
                if not df_axis[name].is_unique:
                    raise KeyError(
                        f"Values in {ax_name}.{name} must be unique. "
                        "Please prepare data to contain unique values."
                    )
                df_axis.reset_index(drop=True, inplace=True)
            else:
                # user specified a non-existent column name
                raise KeyError(
                    f"Annotation name {name}, specified in --{ax_name}-name does not exist."
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
            "annotations": {
                "obs": {
                    "index": self.config["obs_names"],
                    "columns": []
                },
                "var": {
                    "index": self.config["var_names"],
                    "columns": []
                }
            },
            "layout": {"obs": []}
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
                self.schema["annotations"][ax]["columns"].append(ann_schema)

        for layout in self.config['layout']:
            layout_schema = {
                "name": layout,
                "type": "float32",
                "dims": [f"{layout}_0", f"{layout}_1"]
            }
            self.schema["layout"]["obs"].append(layout_schema)

    def _load_data(self, data):
        # as of AnnData 0.6.19, backed mode performs initial load fast, but at the
        # cost of significantly slower access to X data.
        try:
            self.data = anndata.read_h5ad(data)
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
        # var and obs column names must be unique
        if not self.data.obs.columns.is_unique or not self.data.var.columns.is_unique:
            raise KeyError(f"All annotation column names must be unique.")

        self._alias_annotation_names()
        self._validate_data_types()
        self.cell_count = self.data.shape[0]
        self.gene_count = self.data.shape[1]
        self._default_and_validate_layouts()
        self._create_schema()

    @requires_data
    def _default_and_validate_layouts(self):
        """ function:
            a) generate list of default layouts, if not already user specified
            b) validate layouts are legal.  remove/warn on any that are not
            c) cap total list of layouts at global const MAX_LAYOUTS
        """
        layouts = self.config['layout']
        # handle default
        if layouts is None or len(layouts) == 0:
            # load default layouts from the data.
            layouts = [key[2:] for key in self.data.obsm_keys() if type(key) == str and key.startswith("X_")]
            if len(layouts) == 0:
                raise PrepareError(f"Unable to find any precomputed layouts within the dataset.")

        # remove invalid layouts
        valid_layouts = []
        obsm_keys = self.data.obsm_keys()
        for layout in layouts:
            layout_name = f"X_{layout}"
            if layout_name not in obsm_keys:
                warnings.warn(f"Ignoring unknown layout name: {layout}.")
            elif not self._is_valid_layout(self.data.obsm[layout_name]):
                warnings.warn(f"Ignoring layout due to malformed shape or data type: {layout}")
            else:
                valid_layouts.append(layout)

        if len(valid_layouts) == 0:
            raise PrepareError(f"No valid layout data.")

        # cap layouts to MAX_LAYOUTS
        self.config['layout'] = valid_layouts[0:MAX_LAYOUTS]

    @requires_data
    def _is_valid_layout(self, arr):
        """ return True if this layout data is a valid array for front-end presentation:
            * ndarray, with shape (n_obs, >= 2), dtype float/int/uint
            * contains only finite values
        """
        is_valid = type(arr) == np.ndarray and arr.dtype.kind in "fiu"
        is_valid = is_valid and arr.shape[0] == self.data.n_obs and arr.shape[1] >= 2
        is_valid = is_valid and np.all(np.isfinite(arr))
        return is_valid

    @requires_data
    def _validate_data_types(self):
        if sparse.isspmatrix(self.data.X) and not sparse.isspmatrix_csc(self.data.X):
            warnings.warn(
                f"Scanpy data matrix is sparse, but not a CSC (columnar) matrix.  "
                f"Performance may be improved by using CSC."
            )
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

        All embeddings must be individually centered & scaled (isotropically)
        to a [0, 1] range.
        """
        try:
            layout_data = []
            for layout in self.config["layout"]:
                full_embedding = self.data.obsm[f"X_{layout}"]
                embedding = full_embedding[:, :2]

                # scale isotropically
                min = embedding.min(axis=0)
                max = embedding.max(axis=0)
                scale = np.amax(max - min)
                normalized_layout = (embedding - min) / scale

                # translate to center on both axis
                translate = 0.5 - ((max - min) / scale / 2)
                normalized_layout = normalized_layout + translate

                normalized_layout = normalized_layout.astype(dtype=np.float32)
                layout_data.append(pandas.DataFrame(normalized_layout, columns=[f"{layout}_0", f"{layout}_1"]))

        except ValueError as e:
            raise PrepareError(
                f"Layout has not been calculated using {self.config['layout']}, "
                f"please prepare your datafile and relaunch cellxgene") from e

        df = pandas.concat(layout_data, axis=1, copy=False)
        return encode_matrix_fbs(df, col_idx=df.columns, row_idx=None)
