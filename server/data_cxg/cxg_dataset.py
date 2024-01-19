import json
import logging
import os
import threading

import numpy as np
import pandas as pd
import tiledb
from packaging import version
from server_timing import Timing as ServerTiming
from tiledb import TileDBError

from server.common.constants import Axis, XApproximateDistribution
from server.common.errors import ConfigurationError, DatasetAccessError
from server.common.errors import FilterError, JSONEncodingValueError, ExceedsLimitError
from server.common.fbs.matrix import encode_matrix_fbs
from server.common.immutable_kvcache import ImmutableKVCache
from server.common.utils.type_conversion_utils import get_schema_type_hint_from_dtype
from server.common.utils.utils import path_join, jsonify_numpy
from server.compute import diffexp_cxg
from server.data_cxg.cxg_util import pack_selector_from_mask
from server.data_common.data_adaptor import DataAdaptor

class CxgDataset(DataAdaptor):
    # These defaults are overridden by the config variable: server.adaptor.cxg_adaptor.tiledb_cxt
    tiledb_ctx = tiledb.Ctx(
        {
            "sm.tile_cache_size": 8 * 1024**3,
            "py.init_buffer_bytes": 512 * 1024**2,
            "vfs.s3.region": "us-east-1",
        }
    )

    def __init__(self, data_locator, app_config=None, dataset_config=None):
        super().__init__(data_locator, app_config, dataset_config)
        self.lock = threading.Lock()

        self.url = data_locator.uri_or_path
        if self.url[-1] != "/":
            self.url += "/"

        # caching immutable state
        self.lsuri_results = ImmutableKVCache(lambda key: self._lsuri(uri=key, tiledb_ctx=self.tiledb_ctx))
        self.arrays = ImmutableKVCache(lambda key: self._open_array(uri=key, tiledb_ctx=self.tiledb_ctx))
        self.schema = None
        self.X_approximate_distribution = None

        self._validate_and_initialize()

    def cleanup(self):
        """close all the open tiledb arrays"""
        for array in self.arrays.values():
            array.close()
        self.arrays.clear()

    @staticmethod
    def set_tiledb_context(context_params):
        """Set the tiledb context.  This should be set before any instances of CxgDataset are created"""
        try:
            CxgDataset.tiledb_ctx = tiledb.Ctx(context_params)

        except tiledb.libtiledb.TileDBError as e:
            if e.message == "Global context already initialized!":
                if tiledb.default_ctx().config().dict() != CxgDataset.tiledb_ctx.config().dict():
                    raise ConfigurationError("Cannot change tiledb configuration once it is set") from None
            else:
                raise ConfigurationError(f"Invalid tiledb context: {str(e)}") from None

    @staticmethod
    def pre_load_validation(data_locator):
        location = data_locator.uri_or_path
        if not CxgDataset.isvalid(location):
            logging.error(f"cxg matrix is not valid: {location}")
            raise DatasetAccessError("cxg matrix is not valid")

    @staticmethod
    def file_size(data_locator):
        return 0

    @staticmethod
    def open(data_locator, app_config, dataset_config=None):
        return CxgDataset(data_locator, app_config, dataset_config)

    def get_about(self):
        return self.about if self.about else super().get_about()

    def get_title(self):
        return self.title if self.title else super().get_title()

    def get_corpora_props(self):
        return self.corpora_props if self.corpora_props else super().get_corpora_props()

    def get_name(self):
        return "cellxgene cxg adaptor version"

    def get_library_versions(self):
        return dict(tiledb=tiledb.version.version)

    def get_path(self, *urls):
        return path_join(self.url, *urls)

    def _init_sparsity_and_format(self):
        try:
            X = self.open_array("Xr")
            self.is_1d = True
        except Exception:
            X = self.open_array("X")
            self.is_1d = False
        self.is_sparse = X.schema.sparse

    @staticmethod
    def _lsuri(uri, tiledb_ctx):
        def _cleanpath(p):
            if p[-1] == "/":
                return p[:-1]
            else:
                return p

        result = []
        tiledb.ls(uri, lambda path, type: result.append((_cleanpath(path), type)), ctx=tiledb_ctx)
        return result

    def lsuri(self, uri):
        """
        given a URI, do a tiledb.ls but normalizing for all path weirdness:
        * S3 URIs require trailing slash.  file: doesn't care.
        * results on S3 *have* a trailing slash, Posix does not.

        returns list of (absolute paths, type) *without* trailing slash
        in the path.
        """
        if uri[-1] != "/":
            uri += "/"
        return self.lsuri_results[uri]

    @staticmethod
    def isvalid(url):
        """
        Return True if this looks like a valid CXG, False if not.  Just a quick/cheap
        test, not to be fully trusted.
        """
        if tiledb.object_type(url, ctx=CxgDataset.tiledb_ctx) != "group":
            return False
        if tiledb.object_type(path_join(url, "obs"), ctx=CxgDataset.tiledb_ctx) != "array":
            return False
        if tiledb.object_type(path_join(url, "var"), ctx=CxgDataset.tiledb_ctx) != "array":
            return False
        if tiledb.object_type(path_join(url, "X"), ctx=CxgDataset.tiledb_ctx) != "array" and not (
            tiledb.object_type(path_join(url, "Xr"), ctx=CxgDataset.tiledb_ctx) == "array"
            and tiledb.object_type(path_join(url, "Xc"), ctx=CxgDataset.tiledb_ctx) == "array"
        ):
            return False
        if tiledb.object_type(path_join(url, "emb"), ctx=CxgDataset.tiledb_ctx) != "group":
            return False
        return True

    def has_array(self, name):
        a_type = tiledb.object_type(path_join(self.url, name), ctx=self.tiledb_ctx)
        return a_type == "array"

    def _validate_and_initialize(self):
        """
        remember, preload_validation() has already been called, so
        no need to repeat anything it has done.

        Load the CXG "group" metadata and cache instance values.
        Be very aware of multiple versions of the CXG object.

        CXG versions in the wild:
        * version 0, aka "no version" -- can be detected by the lack
          of a cxg_group_metadata array.
        * version 0.1 -- metadata attache to cxg_group_metadata array.
          Same as 0, except it adds group metadata.
        """
        title = None
        about = None
        corpora_props = None
        if self.has_array("cxg_group_metadata"):
            # version >0
            gmd = self.open_array("cxg_group_metadata")
            cxg_version = gmd.meta["cxg_version"]
            if version.parse(cxg_version) >= version.parse("0.1.0"):
                cxg_properties = json.loads(gmd.meta["cxg_properties"])
                title = cxg_properties.get("title", None)
                about = cxg_properties.get("about", None)
            if version.parse(cxg_version) >= version.parse("0.2.0"):
                corpora_props = json.loads(gmd.meta["corpora"]) if "corpora" in gmd.meta else None
        else:
            # version 0
            cxg_version = "0.0"

        if cxg_version not in ["0.0", "0.1", "0.2.0", "0.3.0"]:
            raise DatasetAccessError(f"cxg matrix is not valid: {self.url}")

        self.X_approximate_distribution = self.app_config.dataset_config.X_approximate_distribution

        self.title = title
        self.about = about
        self.cxg_version = cxg_version
        self.corpora_props = corpora_props
        self._init_sparsity_and_format()

    
        annotations = self.get_schema()['annotations']
        self.parameters = {
            'obs_names': annotations['obs']['index'],
            'var_names': annotations['var']['index']
        }

    def data_frame_to_fbs_matrix(self, filter, axis, num_bins=None):
        """
        Retrieves data 'X' and returns in a flatbuffer Matrix.
        :param filter: filter: dictionary with filter params
        :param axis: string obs or var
        :param num_bins: number of bins for lossy integer compression. if None, no compression is performed.
        :return: flatbuffer Matrix

        Caveats:
        * currently only supports access on VAR axis
        * currently only supports filtering on VAR axis
        """
        with ServerTiming.time("where.query"):
            if axis != Axis.VAR:
                raise ValueError("Only VAR dimension access is supported")

            try:
                obs_selector, var_selector = self._filter_to_mask(filter)
            except (KeyError, IndexError, TypeError, AttributeError, DatasetAccessError):
                raise FilterError("Error parsing filter") from None

            if obs_selector is not None:
                raise FilterError("filtering on obs unsupported")

            num_columns = self.get_shape()[1] if var_selector is None else np.count_nonzero(var_selector)
            if self.server_config.exceeds_limit("column_request_max", num_columns):
                raise ExceedsLimitError("Requested dataframe columns exceed column request limit")

            X = self.get_X_array(obs_selector, var_selector)
        with ServerTiming.time("where.encode"):
            col_idx = np.nonzero([] if var_selector is None else var_selector)[0]
            fbs = encode_matrix_fbs(X, col_idx=col_idx, row_idx=None, num_bins=num_bins)

        return fbs

    @staticmethod
    def _open_array(uri, tiledb_ctx):
        return tiledb.open(uri, mode="r", ctx=tiledb_ctx)

    def open_X_array(self, col_wise=False):
        """
        A helper function to open the 1D X array in row or column orientation with
        backwards compatibility for 2D arrays.
        """
        if self.is_1d and col_wise:
            return self.open_array("Xc")
        elif self.is_1d and not col_wise:
            return self.open_array("Xr")
        else:
            # If not 1D, then X array is either dense data or sparse data that has not
            # been remastered. Support for 2D sparse arrays is deprecated.
            return self.open_array("X")

    def open_array(self, name):
        try:
            p = self.get_path(name)
            return self.arrays[p]
        except tiledb.libtiledb.TileDBError:
            raise DatasetAccessError(name) from None

    def get_embedding_array(self, ename, dims=2):
        array = self.open_array(f"emb/{ename}")
        return array[:, 0:dims]

    def diffexp_topN_from_list(self, listA: np.ndarray, listB: np.ndarray, top_n: int = None):
        """
        Compute differential expression - same as diffexp_topN() - but specifying the
        two cell sets as lists of obs indices (postings lists).
        """
        if top_n is None:
            top_n = self.app_config.dataset_config.diffexp__top_n

        result = self.compute_diffexp_ttest(
            listA,
            listB,
            top_n=top_n,
            lfc_cutoff=self.app_config.dataset_config.diffexp__lfc_cutoff,
            selector_lists=True,
        )

        try:
            return jsonify_numpy(result)
        except ValueError:
            raise JSONEncodingValueError("Error encoding differential expression to JSON") from None

    def compute_diffexp_ttest(self, setA, setB, top_n=None, lfc_cutoff=None, selector_lists=False):
        if top_n is None:
            top_n = self.app_config.default_dataset.diffexp__top_n
        if lfc_cutoff is None:
            lfc_cutoff = self.app_config.default_dataset.diffexp__lfc_cutoff
        return diffexp_cxg.diffexp_ttest(
            adaptor=self,
            setA=setA,
            setB=setB,
            top_n=top_n,
            diffexp_lfc_cutoff=lfc_cutoff,
            selector_lists=selector_lists,
        )

    def get_colors(self):
        if self.cxg_version == "0.0":
            return dict()
        meta = self.open_array("cxg_group_metadata").meta
        return json.loads(meta["cxg_category_colors"]) if "cxg_category_colors" in meta else dict()

    def __remap_indices(self, coord_range, coord_mask, coord_data):
        """
        This function maps the indices in coord_data, which could be in the range [0,coord_range), to
        a range that only includes the number of indices encoded in coord_mask.
        coord_range is the maxinum size of the range (e.g. get_shape()[0] or get_shape()[1])
        coord_mask is a mask passed into the get_X_array, of size coord_range
        coord_data are indices representing locations of non-zero values, in the range [0,coord_range).

        For example, say
        coord_mask = [1,0,1,0,0,1]
        coord_data = [2,0,2,2,5]

        The function computes the following:
        indices = [0,2,5]
        ncoord = 3
        maprange = [0,1,2]
        mapindex = [0,0,1,0,0,2]
        coordindices = [1,0,1,1,2]
        """
        if coord_mask is None:
            return coord_range, coord_data

        indices = coord_mask.nonzero()[0]
        ncoord = indices.shape[0]
        maprange = np.arange(ncoord)
        mapindex = np.zeros(indices[-1] + 1, dtype=int)
        mapindex[indices] = maprange
        coordindices = mapindex[coord_data]
        return ncoord, coordindices

    def get_X_array(self, obs_mask=None, var_mask=None):
        obs_items = pack_selector_from_mask(obs_mask)
        var_items = pack_selector_from_mask(var_mask)
        shape = self.get_shape()
        if obs_items is None or var_items is None:
            # If either zero rows or zero columns were selected, return an empty 2d array.
            obs_size = 0 if obs_items is None else shape[0] if obs_mask is None else np.count_nonzero(obs_mask)
            var_size = 0 if var_items is None else shape[1] if var_mask is None else np.count_nonzero(var_mask)
            return np.ndarray((obs_size, var_size))

        if self.is_sparse:
            X = self.open_X_array(col_wise=True)
            if self.is_1d:
                data = X.query(order="U").multi_index[var_items]
            else:
                data = X.query(order="U").multi_index[obs_items, var_items]
            nrows, obsindices = self.__remap_indices(shape[0], obs_mask, data.get("coords", data)["obs"])
            ncols, varindices = self.__remap_indices(shape[1], var_mask, data.get("coords", data)["var"])
            densedata = np.zeros((nrows, ncols), dtype=self.get_X_array_dtype())
            densedata[obsindices, varindices] = data[""]
            return densedata
        else:
            X = self.open_X_array()
            data = X.multi_index[obs_items, var_items][""]
            return data

    def get_X_approximate_distribution(self) -> XApproximateDistribution:
        return self.X_approximate_distribution

    def get_shape(self):
        X = self.open_X_array()
        if self.is_1d:
            Xc = self.open_X_array(col_wise=True)
            return (X.shape[0], Xc.shape[0])
        else:
            return X.shape

    def get_X_array_dtype(self):
        X = self.open_X_array()
        return X.attr(0).dtype

    def query_var_array(self, term_name):
        var = self.open_array("var")
        data = var.query(attrs=[term_name])[:][term_name]
        return data

    def query_obs_array(self, term_name):
        obs = self.open_array("obs")
        schema = self.get_schema()
        try:
            data = obs.query(attrs=[term_name])[:][term_name]
            type_hint = schema.get(term_name, None)
            if (
                type_hint is not None
                and type_hint[term_name]["type"] == "categorical"
                and str(data.dtype).startswith("int")
                and "categories" in type_hint[term_name]
            ):
                data = pd.Categorical.from_codes(data, categories=type_hint[term_name]["categories"])
            elif type_hint is not None and type_hint[term_name]["type"] == "boolean":
                # convert boolean to categorical
                data = pd.Categorical(data)
        except tiledb.libtiledb.TileDBError:
            raise DatasetAccessError("query_obs") from None
        return data

    def get_obs_names(self):
        # get the index from the meta data
        obs = self.open_array("obs")
        meta = json.loads(obs.meta["cxg_schema"])
        index_name = meta["index"]
        return index_name

    def get_obs_index(self):
        obs = self.open_array("obs")
        meta = json.loads(obs.meta["cxg_schema"])
        index_name = meta["index"]
        data = obs.query(attrs=[index_name])[:][index_name]
        return data

    def get_obs_columns(self):
        obs = self.open_array("obs")
        schema = obs.schema
        col_names = [attr.name for attr in schema]
        return pd.Index(col_names)

    def get_obs_keys(self):
        obs = self.open_array("obs")
        schema = obs.schema
        return [attr.name for attr in schema]

    def get_var_keys(self):
        var = self.open_array("var")
        schema = var.schema
        return [attr.name for attr in schema]

    # function to get the embedding
    # this function to iterate through embeddings.
    def get_embedding_names(self):
        with ServerTiming.time("layout.lsuri"):
            pemb = self.get_path("emb")
            embeddings = [os.path.basename(p) for (p, t) in self.lsuri(pemb) if t == "array"]
        if len(embeddings) == 0:
            raise DatasetAccessError("cxg matrix missing embeddings")
        return embeddings

    def _get_schema(self):
        if self.schema:
            return self.schema

        shape = self.get_shape()
        dtype = self.get_X_array_dtype()

        dataframe = {
            "nObs": shape[0],
            "nVar": shape[1],
            # Allow int64 fields to be generated in the schema hint so that we can filter later
            **get_schema_type_hint_from_dtype(dtype=dtype, allow_int64=True),
        }

        annotations = {}
        for ax in ("obs", "var"):
            A = self.open_array(ax)
            schema_hints = json.loads(A.meta["cxg_schema"]) if "cxg_schema" in A.meta else {}
            if not isinstance(schema_hints, dict):
                raise TypeError("Array schema was malformed.")

            cols = []
            for attr in A.schema:
                schema = dict(name=attr.name, writable=False)
                type_hint = schema_hints.get(attr.name, {})
                # type hints take precedence
                if "type" in type_hint:
                    if type_hint["type"] in ["int64", "uint64"] and ax == "obs":
                        # Skip over int64 fields in the obs array when generating schema
                        continue

                    schema["type"] = type_hint["type"]
                    if schema["type"] == "boolean" and ax == "obs":
                        # convert boolean to categorical
                        schema["type"] = "categorical"
                        schema["categories"] = pd.Categorical(
                            self.open_array("obs").query(attrs=[attr.name])[:][attr.name].astype("bool")
                        ).categories.tolist()
                    elif schema["type"] == "categorical" and "categories" in type_hint:
                        schema["categories"] = type_hint["categories"]
                else:
                    schema.update(get_schema_type_hint_from_dtype(dtype=attr.dtype))
                cols.append(schema)

            annotations[ax] = dict(columns=cols)

            if "index" in schema_hints:
                annotations[ax].update({"index": schema_hints["index"]})

        obs_layout = []
        embeddings = self.get_embedding_names()
        for ename in embeddings:
            A = self.open_array(f"emb/{ename}")
            obs_layout.append({"name": ename, "type": "float32", "dims": [f"{ename}_{d}" for d in range(0, A.ndim)]})

        schema = {"dataframe": dataframe, "annotations": annotations, "layout": {"obs": obs_layout}}
        return schema

    def get_schema(self):
        if self.schema is None:
            with self.lock:
                self.schema = self._get_schema()
        return self.schema

    def annotation_to_fbs_matrix(self, axis, fields=None, num_bins=None):
        with ServerTiming.time(f"annotations.{axis}.query"):
            A = self.open_array(str(axis))

            try:
                data = A[:] if not fields else A.query(attrs=fields)[:]

                categorical_dtypes = []
                for c in self.get_schema()["annotations"]["obs"]["columns"]:
                    if c["name"] in fields and c["type"] == "categorical":
                        categories = c.get("categories", None)
                        if categories:
                            categorical_dtypes.append((c["name"], c["categories"]))

                df = pd.DataFrame.from_dict(data)
                for name, categories in categorical_dtypes:
                    if str(df[name].dtype).startswith("int"):
                        df[name] = pd.Categorical.from_codes(df[name], categories=categories)
                    else:
                        df[name] = pd.Categorical(df[name])

                if fields:
                    df = df[fields]
            except TileDBError as e:
                raise KeyError(e) from None

        with ServerTiming.time(f"annotations.{axis}.encode"):
            fbs = encode_matrix_fbs(df, col_idx=df.columns, num_bins=num_bins)

        return fbs
