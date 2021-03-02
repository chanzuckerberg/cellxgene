import json
import logging
import os
import threading

import numpy as np
import pandas as pd
import tiledb
from server_timing import Timing as ServerTiming

import server.compute.diffexp_cxg as diffexp_cxg
from server.common.constants import Axis
from server.common.errors import DatasetAccessError, ConfigurationError
from server.common.immutable_kvcache import ImmutableKVCache
from shared_utils.utils.type_conversion_utils import get_schema_type_hint_from_dtype
from shared_utils.utils.utils import path_join
from server.data_common.data_adaptor import DataAdaptor
from server.data_common.fbs.matrix import encode_matrix_fbs
from server.data_cxg.cxg_util import pack_selector_from_mask


class CxgAdaptor(DataAdaptor):
    # TODO:  The tiledb context parameters should be a configuration option
    tiledb_ctx = tiledb.Ctx(
        {"sm.tile_cache_size": 8 * 1024 * 1024 * 1024, "sm.num_reader_threads": 32, "vfs.s3.region": "us-east-1"}
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

        self._validate_and_initialize()

    def cleanup(self):
        """close all the open tiledb arrays"""
        for array in self.arrays.values():
            array.close()
        self.arrays.clear()

    @staticmethod
    def set_tiledb_context(context_params):
        """Set the tiledb context.  This should be set before any instances of CxgAdaptor are created"""
        try:
            CxgAdaptor.tiledb_ctx = tiledb.Ctx(context_params)
            tiledb.default_ctx(context_params)

        except tiledb.libtiledb.TileDBError as e:
            if e.message == "Global context already initialized!":
                if tiledb.default_ctx().config().dict() != CxgAdaptor.tiledb_ctx.config().dict():
                    raise ConfigurationError("Cannot change tiledb configuration once it is set")
            else:
                raise ConfigurationError(f"Invalid tiledb context: {str(e)}")

    @staticmethod
    def pre_load_validation(data_locator):
        location = data_locator.uri_or_path
        if not CxgAdaptor.isvalid(location):
            logging.error(f"cxg matrix is not valid: {location}")
            raise DatasetAccessError("cxg matrix is not valid")

    @staticmethod
    def file_size(data_locator):
        return 0

    @staticmethod
    def open(data_locator, app_config, dataset_config=None):
        return CxgAdaptor(data_locator, app_config, dataset_config)

    def get_about(self):
        return self.about if self.about else super().get_about()

    def get_title(self):
        return self.title if self.title else super().get_title()

    def get_corpora_props(self):
        return self.corpora_props if self.corpora_props else super().get_corpora_props()

    def get_name(self):
        return "cellxgene cxg adaptor version"

    def get_library_versions(self):
        return dict(tiledb=tiledb.__version__)

    def get_path(self, *urls):
        return path_join(self.url, *urls)

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
        if not tiledb.object_type(url, ctx=CxgAdaptor.tiledb_ctx) == "group":
            return False
        if not tiledb.object_type(path_join(url, "obs"), ctx=CxgAdaptor.tiledb_ctx) == "array":
            return False
        if not tiledb.object_type(path_join(url, "var"), ctx=CxgAdaptor.tiledb_ctx) == "array":
            return False
        if not tiledb.object_type(path_join(url, "X"), ctx=CxgAdaptor.tiledb_ctx) == "array":
            return False
        if not tiledb.object_type(path_join(url, "emb"), ctx=CxgAdaptor.tiledb_ctx) == "group":
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
            # version 0.1 used a malformed/shorthand semver string.
            if cxg_version == "0.1" or cxg_version == "0.2.0":
                cxg_properties = json.loads(gmd.meta["cxg_properties"])
                title = cxg_properties.get("title", None)
                about = cxg_properties.get("about", None)
            if cxg_version == "0.2.0":
                corpora_props = json.loads(gmd.meta["corpora"]) if "corpora" in gmd.meta else None
        else:
            # version 0
            cxg_version = "0.0"

        if cxg_version not in ["0.0", "0.1", "0.2.0"]:
            raise DatasetAccessError(f"cxg matrix is not valid: {self.url}")

        self.title = title
        self.about = about
        self.cxg_version = cxg_version
        self.corpora_props = corpora_props

    @staticmethod
    def _open_array(uri, tiledb_ctx):
        with tiledb.Array(uri, mode="r", ctx=tiledb_ctx) as array:
            if array.schema.sparse:
                return tiledb.SparseArray(uri, mode="r", ctx=tiledb_ctx)
            else:
                return tiledb.DenseArray(uri, mode="r", ctx=tiledb_ctx)

    def open_array(self, name):
        try:
            p = self.get_path(name)
            return self.arrays[p]
        except tiledb.libtiledb.TileDBError:
            raise DatasetAccessError(name)

    def get_embedding_array(self, ename, dims=2):
        array = self.open_array(f"emb/{ename}")
        return array[:, 0:dims]

    def compute_embedding(self, method, filter):
        raise NotImplementedError("CXG does not yet support re-embedding")

    def compute_diffexp_ttest(self, maskA, maskB, top_n=None, lfc_cutoff=None):
        if top_n is None:
            top_n = self.dataset_config.diffexp__top_n
        if lfc_cutoff is None:
            lfc_cutoff = self.dataset_config.diffexp__lfc_cutoff
        return diffexp_cxg.diffexp_ttest(self, maskA, maskB, top_n, lfc_cutoff)

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

        indices = np.where(coord_mask)[0]
        ncoord = indices.shape[0]
        maprange = np.arange(ncoord)
        mapindex = np.zeros(indices[-1] + 1, dtype=int)
        mapindex[indices] = maprange
        coordindices = mapindex[coord_data]
        return ncoord, coordindices

    def get_X_array(self, obs_mask=None, var_mask=None):
        obs_items = pack_selector_from_mask(obs_mask)
        var_items = pack_selector_from_mask(var_mask)
        if obs_items is None or var_items is None:
            # If either zero rows or zero columns were selected, return an empty 2d array.
            shape = self.get_shape()
            obs_size = 0 if obs_items is None else shape[0] if obs_mask is None else np.count_nonzero(obs_mask)
            var_size = 0 if var_items is None else shape[1] if var_mask is None else np.count_nonzero(var_mask)
            return np.ndarray((obs_size, var_size))

        X = self.open_array("X")

        if X.schema.sparse:
            if obs_items == slice(None) and var_items == slice(None):
                data = X[:, :]
            else:
                data = X.multi_index[obs_items, var_items]

            nrows, obsindices = self.__remap_indices(X.shape[0], obs_mask, data.get("coords", data)["obs"])
            ncols, varindices = self.__remap_indices(X.shape[1], var_mask, data.get("coords", data)["var"])
            densedata = np.zeros((nrows, ncols), dtype=self.get_X_array_dtype())
            densedata[obsindices, varindices] = data[""]
            if self.has_array("X_col_shift"):
                X_col_shift = self.open_array("X_col_shift")
                if var_items == slice(None):
                    densedata += X_col_shift[:]
                else:
                    densedata += X_col_shift.multi_index[var_items][""]

            return densedata

        else:
            if obs_items == slice(None) and var_items == slice(None):
                data = X[:, :]
            else:
                data = X.multi_index[obs_items, var_items][""]
            return data

    def get_shape(self):
        X = self.open_array("X")
        return X.shape

    def get_X_array_dtype(self):
        X = self.open_array("X")
        return X.dtype

    def query_var_array(self, term_name):
        var = self.open_array("var")
        data = var.query(attrs=[term_name])[:][term_name]
        return data

    def query_obs_array(self, term_name):
        var = self.open_array("obs")
        try:
            data = var.query(attrs=[term_name])[:][term_name]
        except tiledb.libtiledb.TileDBError:
            raise DatasetAccessError("query_obs")
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

        dataframe = {"nObs": shape[0], "nVar": shape[1], "type": dtype.name}

        annotations = {}
        for ax in ("obs", "var"):
            A = self.open_array(ax)
            schema_hints = json.loads(A.meta["cxg_schema"]) if "cxg_schema" in A.meta else {}
            if type(schema_hints) is not dict:
                raise TypeError("Array schema was malformed.")

            cols = []
            for attr in A.schema:
                schema = dict(name=attr.name, writable=False)
                type_hint = schema_hints.get(attr.name, {})
                # type hints take precedence
                if "type" in type_hint:
                    schema["type"] = type_hint["type"]
                    if schema["type"] == "categorical" and "categories" in type_hint:
                        schema["categories"] = type_hint["categories"]
                else:
                    schema.update(get_schema_type_hint_from_dtype(attr.dtype))
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

    def _annotations_field_split(self, axis, fields, A, labels):
        """
        fields: requested fields, may be None (all)
        labels: writable user annotations dataframe, if any

        Remove redundant fields, raise KeyError on non-existant fields,
        and split into three lists:
            fields_to_fetch_from_cxg
            fields_to_fetch_from_labels
            fields_to_return

        if we have to return from labels, the fetch fields will contain the index
        to join on, which may not be in fields_to_return
        """
        need_labels = axis == Axis.OBS and labels is not None and not labels.empty
        index_key = self.get_obs_names() if need_labels else None

        if not fields:
            return (None, None, None, index_key)

        cxg_keys = frozenset([a.name for a in A.schema])
        user_anno_keys = frozenset(labels.columns.tolist()) if need_labels else frozenset()
        return_keys = frozenset(fields)

        label_join_index = frozenset([index_key]) if need_labels and (return_keys & user_anno_keys) else frozenset()

        unknown_fields = return_keys - (cxg_keys | user_anno_keys)
        if unknown_fields:
            raise KeyError("_".join(unknown_fields))

        return (
            list((return_keys & cxg_keys) | label_join_index),
            list(return_keys & user_anno_keys),
            list(return_keys),
            index_key,
        )

    def annotation_to_fbs_matrix(self, axis, fields=None, labels=None):
        with ServerTiming.time(f"annotations.{axis}.query"):
            A = self.open_array(str(axis))

            # may raise if fields contains unknown key
            cxg_fields, anno_fields, return_fields, index_field = self._annotations_field_split(axis, fields, A, labels)

            if cxg_fields is None:
                data = A[:]
            elif cxg_fields:
                data = A.query(attrs=cxg_fields)[:]
            else:
                data = {}

            df = pd.DataFrame.from_dict(data)

            if axis == Axis.OBS and labels is not None and not labels.empty:
                if anno_fields is None:
                    assert index_field
                    df = df.join(labels, index_field)
                elif anno_fields:
                    assert index_field
                    df = df.join(labels[anno_fields], index_field)

            if return_fields:
                df = df[return_fields]

        with ServerTiming.time(f"annotations.{axis}.encode"):
            fbs = encode_matrix_fbs(df, col_idx=df.columns)

        return fbs
