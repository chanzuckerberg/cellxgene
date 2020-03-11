import warnings

import numpy as np
import pandas as pd
from pandas.core.dtypes.dtypes import CategoricalDtype
import anndata
from scipy import sparse
from packaging import version
from datetime import datetime
from server_timing import Timing as ServerTiming

from server.data_common.data_adaptor import DataAdaptor
from server.data_common.fbs.matrix import encode_matrix_fbs
from server.common.utils import series_to_schema
from server.common.constants import Axis, MAX_LAYOUTS
from server.common.errors import PrepareError, DatasetAccessError, FilterError
from server.compute.scanpy import scanpy_umap

anndata_version = version.parse(str(anndata.__version__)).release


def anndata_version_is_pre_070():
    major = anndata_version[0]
    minor = anndata_version[1] if len(anndata_version) > 1 else 0
    return major == 0 and minor < 7


class AnndataAdaptor(DataAdaptor):
    def __init__(self, data_locator, config=None):
        super().__init__(config)
        self.data = None
        self.data_locator = data_locator
        self._load_data(data_locator)
        self._validate_and_initialize()

    def cleanup(self):
        pass

    @staticmethod
    def pre_load_validation(data_locator):
        if data_locator.islocal():
            # if data locator is local, apply file system conventions and other "cheap"
            # validation checks.  If a URI, defer until we actually fetch the data and
            # try to read it.  Many of these tests don't make sense for URIs (eg, extension-
            # based typing).
            if not data_locator.exists():
                raise DatasetAccessError(f"{data_locator.uri_or_path} does not exist")
            if not data_locator.isfile():
                raise DatasetAccessError(f"{data_locator.uri_or_path} is not a file")

    @staticmethod
    def file_size(data_locator):
        return data_locator.size() if data_locator.islocal() else 0

    @staticmethod
    def open(data_locator, config):
        return AnndataAdaptor(data_locator, config)

    def get_location(self):
        return self.data_locator.uri_or_path

    def get_data_locator(self):
        return self.data_locator

    def get_name(self):
        return "cellxgene anndata adaptor version"

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
                self.parameters[config_name] = name
            else:
                # user specified a non-existent column name
                raise KeyError(f"Annotation name {name}, specified in --{ax_name}-name does not exist.")

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
                ann_schema.update(series_to_schema(curr_axis[ann]))
                self.schema["annotations"][ax]["columns"].append(ann_schema)

        for layout in self.get_embedding_names():
            layout_schema = {"name": layout, "type": "float32", "dims": [f"{layout}_0", f"{layout}_1"]}
            self.schema["layout"]["obs"].append(layout_schema)

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
                backed = "r" if self.config.anndata_backed else None
                self.data = anndata.read_h5ad(lh, backed=backed)

        except ValueError:
            raise DatasetAccessError(
                "File must be in the .h5ad format. Please read "
                "https://github.com/theislab/scanpy_usage/blob/master/170505_seurat/info_h5ad.md to "
                "learn more about this format. You may be able to convert your file into this format "
                "using `cellxgene prepare`, please run `cellxgene prepare --help` for more "
                "information."
            )
        except MemoryError:
            raise DatasetAccessError("Out of memory - file is too large for available memory.")
        except Exception as e:
            raise DatasetAccessError(
                f"{e} - file not found or is inaccessible.  File must be an .h5ad object.  "
                f"Please check your input and try again."
            )

    def _validate_and_initialize(self):
        if anndata_version_is_pre_070() and self.config.anndata_backed:
            warnings.warn(
                f"Use of --backed mode with anndata versions older than 0.7 will have serious "
                "performance issues. Please update to at least anndata 0.7 or later."
            )

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
        if (n_values > 1e8 and self.config.anndata_backed is True) or (n_values > 5e8):
            self.parameters.update({"diffexp_may_be_slow": True})

    def _is_valid_layout(self, arr):
        """ return True if this layout data is a valid array for front-end presentation:
            * ndarray, with shape (n_obs, >= 2), dtype float/int/uint
            * contains only finite values
        """
        is_valid = type(arr) == np.ndarray and arr.dtype.kind in "fiu"
        is_valid = is_valid and arr.shape[0] == self.data.n_obs and arr.shape[1] >= 2
        is_valid = is_valid and np.all(np.isfinite(arr))
        return is_valid

    def _validate_data_types(self):
        # The backed API does not support interrogation of the underlying sparsity or sparse matrix type
        # Fake it by asking for a small subarray and testing it.   NOTE: if the user has ignored our
        # anndata <= 0.7 warning, opted for the --backed option, and specified a large, sparse dataset,
        # this "small" indexing request will load the entire X array. This is due to a bug in anndata<=0.7
        # which will load the entire X matrix to fullfill any slicing request if X is sparse.  See
        # user warning in _load_data().
        X0 = self.data.X[0, 0:1]
        if sparse.isspmatrix(X0) and not sparse.isspmatrix_csc(X0):
            warnings.warn(
                f"Anndata data matrix is sparse, but not a CSC (columnar) matrix.  "
                f"Performance may be improved by using CSC."
            )
        if self.data.X.dtype != "float32":
            warnings.warn(
                f"Anndata data matrix is in {self.data.X.dtype} format not float32. " f"Precision may be truncated."
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
                        f"Anndata annotation {ax}:{ann} is in unsupported format: {datatype}. "
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
        """
        Return pre-computed embeddings.

        function:
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

    def get_embedding_array(self, ename, dims=2):
        full_embedding = self.data.obsm[f"X_{ename}"]
        return full_embedding[:, 0:dims]

    def compute_embedding(self, method, obsFilter):
        if Axis.VAR in obsFilter:
            raise FilterError("Observation filters may not contain variable conditions")
        if method != "umap":
            raise NotImplementedError(f"re-embedding method {method} is not available.")
        try:
            shape = self.get_shape()
            obs_mask = self._axis_filter_to_mask(Axis.OBS, obsFilter["obs"], shape[0])
        except (KeyError, IndexError) as e:
            raise FilterError(f"Error parsing filter: {e}") from e

        with ServerTiming.time("layout.compute"):
            X_umap = scanpy_umap(self.data, obs_mask)
            normalized_layout = DataAdaptor.normalize_embedding(X_umap)

        # Server picks reemedding name, which must not collide with any other
        # embedding name generated by this backed.
        name = f"reembed:{method}_{datetime.now().isoformat(timespec='milliseconds')}"
        dims = [f"{name}_0", f"{name}_1"]
        df = pd.DataFrame(normalized_layout, columns=dims)
        fbs = encode_matrix_fbs(df, col_idx=df.columns, row_idx=None)
        schema = {"name": name, "type": "float32", "dims": dims}
        return (schema, fbs)

    def get_X_array(self, obs_mask=None, var_mask=None):
        if obs_mask is None:
            obs_mask = slice(None)
        if var_mask is None:
            var_mask = slice(None)
        X = self.data.X[obs_mask, var_mask]
        return X

    def get_shape(self):
        return self.data.shape

    def query_var_array(self, term_name):
        return getattr(self.data.var, term_name)

    def query_obs_array(self, term_name):
        return getattr(self.data.obs, term_name)

    def get_obs_index(self):
        name = getattr(self.config, "obs_names")
        if name is None:
            return self.original_obs_index
        else:
            return self.data.obs[name]

    def get_obs_columns(self):
        return self.data.obs.columns
