import warnings

import numpy as np
from pandas.core.dtypes.dtypes import CategoricalDtype
import anndata
from scipy import sparse

from server.data_common.driver import CXGDriver
from server.common.constants import Axis, MAX_LAYOUTS
from server.common.errors import (
    PrepareError,
    ScanpyFileError,
)
from server.data_common.utils import requires_data
from server.data_common.fbs.matrix import encode_matrix_fbs
import server.data_scanpy.matrix_proxy  # noqa: F401

from server.common.data_locator import DataLocator


class ScanpyEngine(CXGDriver):
    def __init__(self, data_locator=None, config=None):
        super().__init__(config)
        self.data = None
        self.data_locator = None
        self.update(data_locator, config)

    def update(self, data_locator=None, config=None):
        super().__init__(config)
        if data_locator:
            self.data_locator = data_locator
            self._load_data(data_locator)
            self._validate_and_initialize()

    @staticmethod
    def pre_check(location):
        data_locator = DataLocator(location)
        if data_locator.islocal():
            # if data locator is local, apply file system conventions and other "cheap"
            # validation checks.  If a URI, defer until we actually fetch the data and
            # try to read it.  Many of these tests don't make sense for URIs (eg, extension-
            # based typing).
            if not data_locator.exists():
                raise RuntimeError(f"{location} does not exist")
            if not data_locator.isfile():
                raise RuntimeError(f"{location} is not a file")

    @staticmethod
    def file_size(location):
        data_locator = DataLocator(location)
        return data_locator.size() if data_locator.islocal() else 0

    @staticmethod
    def open(location, args):
        data_locator = DataLocator(location)
        return ScanpyEngine(data_locator, args)

    def get_name(self):
        return "cellxgene Scanpy engine version "

    def get_library_versions(self):
        return dict(anndata=str(anndata.__version__))

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
        self.original_obs_index = self.data.obs.index

        for (ax_name, config_name) in ((Axis.OBS, "obs_names"), (Axis.VAR, "var_names")):
            name = getattr(self.config, config_name)
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
                self.parameters[config_name] = name
                # reset index to simple range; alias name to point at the
                # previously specified index.
                df_axis.rename_axis(name, inplace=True)
                df_axis.reset_index(inplace=True)
            elif name in df_axis.columns:
                # User has specified alternative column for unique names, and it exists
                if not df_axis[name].is_unique:
                    raise KeyError(
                        f"Values in {ax_name}.{name} must be unique. " "Please prepare data to contain unique values."
                    )
                df_axis.reset_index(drop=True, inplace=True)
            else:
                # user specified a non-existent column name
                raise KeyError(f"Annotation name {name}, specified in --{ax_name}-name does not exist.")

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

    @staticmethod
    def _get_col_type(col):
        dtype = col.dtype
        data_kind = dtype.kind
        schema = {}

        if ScanpyEngine._can_cast_to_float32(col):
            schema["type"] = "float32"
        elif ScanpyEngine._can_cast_to_int32(col):
            schema["type"] = "int32"
        elif dtype == np.bool_:
            schema["type"] = "boolean"
        elif data_kind == "O" and dtype == "object":
            schema["type"] = "string"
        elif data_kind == "O" and dtype == "category":
            schema["type"] = "categorical"
            schema["categories"] = dtype.categories.tolist()
        else:
            raise TypeError(f"Annotations of type {dtype} are unsupported by cellxgene.")
        return schema

    @requires_data
    def _create_schema(self):
        self.schema = {
            "dataframe": {"nObs": self.cell_count, "nVar": self.gene_count, "type": str(self.data.X.dtype)},
            "annotations": {
                "obs": {"index": self.parameters.get("obs_names"), "columns": []},
                "var": {"index": self.parameters.get("var_names"), "columns": []},
            },
            "layout": {"obs": []},
        }
        for ax in Axis:
            curr_axis = getattr(self.data, str(ax))
            for ann in curr_axis:
                ann_schema = {"name": ann, "writable": False}
                ann_schema.update(self._get_col_type(curr_axis[ann]))
                self.schema["annotations"][ax]["columns"].append(ann_schema)

        for layout in self.get_embedding_names():
            layout_schema = {"name": layout, "type": "float32", "dims": [f"{layout}_0", f"{layout}_1"]}
            self.schema["layout"]["obs"].append(layout_schema)

    @requires_data
    def get_schema(self):
        return self.schema

    def _load_data(self, data_locator):
        # as of AnnData 0.6.19, backed mode performs initial load fast, but at the
        # cost of significantly slower access to X data.
        try:
            # there is no guarantee data_locator indicates a local file.  The AnnData
            # API will only consume local file objects.  If we get a non-local object,
            # make a copy in tmp, and delete it after we load into memory.
            with data_locator.local_handle() as lh:
                # as of AnnData 0.6.19, backed mode performs initial load fast, but at the
                # cost of significantly slower access to X data.
                backed = "r" if self.config.scanpy_backed else None
                self.data = anndata.read_h5ad(lh, backed=backed)

        except ValueError:
            raise ScanpyFileError(
                "File must be in the .h5ad format. Please read "
                "https://github.com/theislab/scanpy_usage/blob/master/170505_seurat/info_h5ad.md to "
                "learn more about this format. You may be able to convert your file into this format "
                "using `cellxgene prepare`, please run `cellxgene prepare --help` for more "
                "information."
            )
        except MemoryError:
            raise ScanpyFileError("Out of memory - file is too large for available memory.")
        except Exception as e:
            raise ScanpyFileError(
                f"{e} - file not found or is inaccessible.  File must be an .h5ad object.  "
                f"Please check your input and try again."
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
        self._create_schema()

        # heuristic
        n_values = self.data.shape[0] * self.data.shape[1]
        if (n_values > 1e8 and self.config.backed is True) or (n_values > 5e8):
            self.parameters.update({"diffexp_may_be_slow": True})

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
                f"Scanpy data matrix is in {self.data.X.dtype} format not float32. " f"Precision may be truncated."
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
                    if category_num > 500 and category_num > self.config.max_category_items:
                        warnings.warn(
                            f"{str(ax).title()} annotation '{ann}' has {category_num} categories, this may be "
                            f"cumbersome or slow to display. We recommend setting the "
                            f"--max-category-items option to 500, this will hide categorical "
                            f"annotations with more than 500 categories in the UI"
                        )

    @requires_data
    def validate_label_data(self, labels):
        """
        labels is None if disabled, empty if enabled by no data
        """
        if labels is None or labels.empty:
            return

        # all lables must have a name, which must be unique and not used in obs column names
        if not labels.columns.is_unique:
            raise KeyError(f"All column names specified in user annotations must be unique.")

        # the label index must be unique, and must have same values the anndata obs index
        if not labels.index.is_unique:
            raise KeyError(f"All row index values specified in user annotations must be unique.")

        if not labels.index.equals(self.original_obs_index):
            raise KeyError(
                "Label file row index does not match H5AD file index. "
                "Please ensure that column zero (0) in the label file contain the same "
                "index values as the H5AD file."
            )

        duplicate_columns = list(set(labels.columns) & set(self.data.obs.columns))
        if len(duplicate_columns) > 0:
            raise KeyError(
                f"Labels file may not contain column names which overlap " f"with h5ad obs columns {duplicate_columns}"
            )

        # labels must have same count as obs annotations
        if labels.shape[0] != self.data.obs.shape[0]:
            raise ValueError("Labels file must have same number of rows as h5ad file.")

    @requires_data
    def annotation_to_fbs_matrix(self, axis, fields=None, labels=None):
        if axis == Axis.OBS:
            if labels is not None and not labels.empty:
                df = self.data.obs.join(labels, self.parameters.get("obs_names"))
            else:
                df = self.data.obs
        else:
            df = self.data.var

        if fields is not None and len(fields) > 0:
            df = df[fields]
        return encode_matrix_fbs(df, col_idx=df.columns)

    def get_embedding_names(self):
        """ function:
            a) generate list of default layouts
            b) validate layouts are legal.  remove/warn on any that are not
            c) cap total list of layouts at global const MAX_LAYOUTS
        """
        # load default layouts from the data.
        layouts = self.config.layout

        if layouts is None or len(layouts) == 0:
            layouts = [key[2:] for key in self.data.obsm_keys() if type(key) == str and key.startswith("X_")]

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
        return layouts[0:MAX_LAYOUTS]

    def get_embedding_array(self, ename, items=None):
        if items is None:
            items = slice(None)
        full_embedding = self.data.obsm[f"X_{ename}"]
        return full_embedding[items]

    def get_X_array(self, obs_items=None, var_items=None):
        if obs_items is None:
            obs_items = slice(None)
        if var_items is None:
            var_items = slice(None)
        X = self.data.X[obs_items, var_items]
        return X

    def get_X_array_shape(self):
        return self.data.shape

    def query_var_array(self, term_name):
        return getattr(self.data.var, term_name)

    def query_obs_array(self, term_name):
        return getattr(self.data.obs, term_name)
