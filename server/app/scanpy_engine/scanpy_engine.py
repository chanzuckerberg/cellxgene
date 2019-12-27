import warnings
import copy
import threading
from datetime import datetime
import os.path
from hashlib import blake2b
import base64

import numpy as np
import pandas
from pandas.core.dtypes.dtypes import CategoricalDtype
import anndata
from scipy import sparse

from server import __version__ as cellxgene_version
from server.app.driver.driver import CXGDriver
from server.app.util.constants import Axis, DEFAULT_TOP_N, MAX_LAYOUTS
from server.app.util.errors import (
    FilterError,
    JSONEncodingValueError,
    PrepareError,
    ScanpyFileError,
    DisabledFeatureError,
)
from server.app.util.utils import jsonify_scanpy, requires_data
from server.app.scanpy_engine.diffexp import diffexp_ttest
from server.app.util.fbs.matrix import encode_matrix_fbs, decode_matrix_fbs
from server.app.scanpy_engine.labels import read_labels, write_labels
import server.app.scanpy_engine.matrix_proxy  # noqa: F401
from server.app.util.matrix_proxy import MatrixProxy


def has_method(o, name):
    """ return True if `o` has callable method `name` """
    op = getattr(o, name, None)
    return op is not None and callable(op)


class ScanpyEngine(CXGDriver):
    def __init__(self, data_locator=None, args={}):
        super().__init__(data_locator, args)
        # lock used to protect label file write ops
        self.label_lock = threading.RLock()
        if self.data:
            self._validate_and_initialize()

    def update(self, data_locator=None, args={}):
        super().__init__(data_locator, args)
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
            "annotations": False,
            "annotations_file": None,
            "annotations_output_dir": None,
            "backed": False,
            "disable_diffexp": False,
            "diffexp_may_be_slow": False
        }

    def get_config_parameters(self, uid=None, collection=None):
        params = {
            "max-category-items": self.config["max_category_items"],
            "disable-diffexp": self.config["disable_diffexp"],
            "diffexp-may-be-slow": self.config["diffexp_may_be_slow"],
            "annotations": self.config["annotations"]
        }
        if self.config["annotations"]:
            if uid is not None:
                params.update({
                    "annotations-user-data-idhash": self.get_userdata_idhash(uid)
                })
            if self.config['annotations_file'] is not None:
                # user has hard-wired the name of the annotation data collection
                fname = os.path.basename(self.config['annotations_file'])
                collection_fname = os.path.splitext(fname)[0]
                params.update({
                    'annotations-data-collection-is-read-only': True,
                    'annotations-data-collection-name': collection_fname
                })
            elif collection is not None:
                params.update({
                    'annotations-data-collection-is-read-only': False,
                    'annotations-data-collection-name': collection
                })
        return params

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
            raise TypeError(
                f"Annotations of type {dtype} are unsupported by cellxgene."
            )
        return schema

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
                ann_schema = {"name": ann, "writable": False}
                ann_schema.update(self._get_col_type(curr_axis[ann]))
                self.schema["annotations"][ax]["columns"].append(ann_schema)

        for layout in self.config['layout']:
            layout_schema = {
                "name": layout,
                "type": "float32",
                "dims": [f"{layout}_0", f"{layout}_1"]
            }
            self.schema["layout"]["obs"].append(layout_schema)

    @requires_data
    def get_schema(self, uid=None, collection=None):
        schema = self.schema  # base schema
        # add label obs annotations as needed
        labels = read_labels(self.get_anno_fname(uid, collection))
        if labels is not None and not labels.empty:
            schema = copy.deepcopy(schema)
            for col in labels.columns:
                col_schema = {
                    "name": col,
                    "writable": True,
                }
                col_schema.update(self._get_col_type(labels[col]))
                schema["annotations"]["obs"]["columns"].append(col_schema)
        return schema

    def get_userdata_idhash(self, uid):
        """
        Return a short hash that weakly identifies the user and dataset.
        Used to create safe annotations output file names.
        """
        id = (uid + self.data_locator.abspath()).encode()
        idhash = base64.b32encode(blake2b(id, digest_size=5).digest()).decode('utf-8')
        return idhash

    def get_anno_fname(self, uid=None, collection=None):
        """ return the current annotation file name """
        if not self.config["annotations"]:
            return None

        if self.config["annotations_file"] is not None:
            return self.config["annotations_file"]

        # we need to generate a file name, which we can only do if we have a UID and collection name
        if uid is None or collection is None:
            return None
        idhash = self.get_userdata_idhash(uid)
        return os.path.join(self.get_anno_output_dir(), f"{collection}-{idhash}.csv")

    def get_anno_output_dir(self):
        """ return the current annotation output directory """
        if not self.config["annotations"]:
            return None

        if self.config['annotations_output_dir']:
            return self.config['annotations_output_dir']

        if self.config['annotations_file']:
            return os.path.dirname(os.path.abspath(self.config['annotations_file']))

        return os.getcwd()

    def get_anno_backup_dir(self, uid, collection=None):
        """ return the current annotation backup directory """
        if not self.config["annotations"]:
            return None

        fname = self.get_anno_fname(uid, collection)
        root, ext = os.path.splitext(fname)
        return f"{root}-backups"

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
                backed = 'r' if self.config['backed'] else None
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
        self._default_and_validate_layouts()
        self._create_schema()

        # if the user has specified a fixed label file, go ahead and validate it
        # so that we can remove errors early in the process.
        if self.config["annotations_file"]:
            self._validate_label_data(read_labels(self.get_anno_fname()))

        # heuristic
        n_values = self.data.shape[0] * self.data.shape[1]
        if (n_values > 1e8 and self.config['backed'] is True) or (n_values > 5e8):
            self.config.update({"diffexp_may_be_slow": True})

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

    @requires_data
    def _validate_label_data(self, labels):
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
            raise KeyError("Label file row index does not match H5AD file index. "
                           "Please ensure that column zero (0) in the label file contain the same "
                           "index values as the H5AD file.")

        duplicate_columns = list(set(labels.columns) & set(self.data.obs.columns))
        if len(duplicate_columns) > 0:
            raise KeyError(f"Labels file may not contain column names which overlap "
                           f"with h5ad obs columns {duplicate_columns}")

        # labels must have same count as obs annotations
        if labels.shape[0] != self.data.obs.shape[0]:
            raise ValueError("Labels file must have same number of rows as h5ad file.")

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
    def annotation_to_fbs_matrix(self, axis, fields=None, uid=None, collection=None):
        if axis == Axis.OBS:
            if self.config["annotations"]:
                try:
                    labels = read_labels(self.get_anno_fname(uid, collection))
                except Exception as e:
                    raise ScanpyFileError(
                        f"Error while loading label file: {e}, File must be in the .csv format, please check "
                        f"your input and try again."
                    )
            else:
                labels = None

            if labels is not None and not labels.empty:
                df = self.data.obs.join(labels, self.config['obs_names'])
            else:
                df = self.data.obs
        else:
            df = self.data.var
        if fields is not None and len(fields) > 0:
            df = df[fields]
        return encode_matrix_fbs(df, col_idx=df.columns)

    @requires_data
    def annotation_put_fbs(self, axis, fbs, uid=None, collection=None):
        if not self.config["annotations"]:
            raise DisabledFeatureError("Writable annotations are not enabled")

        fname = self.get_anno_fname(uid, collection)
        if not fname:
            raise ScanpyFileError("Writable annotations - unable to determine file name for annotations")

        if axis != Axis.OBS:
            raise ValueError("Only OBS dimension access is supported")

        new_label_df = decode_matrix_fbs(fbs)
        if not new_label_df.empty:
            new_label_df.index = self.original_obs_index
        self._validate_label_data(new_label_df)  # paranoia

        # if any of the new column labels overlap with our existing labels, raise error
        duplicate_columns = list(set(new_label_df.columns) & set(self.data.obs.columns))
        if not new_label_df.columns.is_unique or len(duplicate_columns) > 0:
            raise KeyError(f"Labels file may not contain column names which overlap "
                           f"with h5ad obs columns {duplicate_columns}")

        # update our internal state and save it.  Multi-threading often enabled,
        # so treat this as a critical section.
        with self.label_lock:
            lastmod = self.data_locator.lastmodtime()
            lastmodstr = "'unknown'" if lastmod is None else lastmod.isoformat(timespec="seconds")
            header = f"# Annotations generated on {datetime.now().isoformat(timespec='seconds')} " \
                     f"using cellxgene version {cellxgene_version}\n" \
                     f"# Input data file was {self.data_locator.uri_or_path}, " \
                     f"which was last modified on {lastmodstr}\n"
            write_labels(fname, new_label_df, header, backup_dir=self.get_anno_backup_dir(uid, collection))

        return jsonify_scanpy({"status": "OK"})

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
        X = MatrixProxy.create(self.data.X if var_selector is None
                               else self.data.X[:, var_selector])
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
