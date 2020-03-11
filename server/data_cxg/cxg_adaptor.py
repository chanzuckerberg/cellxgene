import os
import json
from server.common.utils import dtype_to_schema
from server.common.errors import DatasetAccessError
from server.common.utils import path_join
from server.common.constants import Axis
from server.data_common.data_adaptor import DataAdaptor
from server.data_common.fbs.matrix import encode_matrix_fbs
import tiledb
import numpy as np
import pandas as pd
from server_timing import Timing as ServerTiming
import threading


class CxgAdaptor(DataAdaptor):

    # TODO:  The tiledb context parameters should be a configuration option
    tiledb_ctx = tiledb.Ctx({"sm.tile_cache_size": 8 * 1024 * 1024 * 1024, "sm.num_reader_threads": 32})

    def __init__(self, data_locator, config=None):
        super().__init__(config)
        self.arrays = {}
        self.lock = threading.Lock()

        self.data_locator = data_locator
        self.url = data_locator.uri_or_path
        if self.url[-1] != "/":
            self.url += "/"

        self._validate_and_initialize()

    def cleanup(self):
        """close all the open tiledb arrays"""
        for array in self.arrays.values():
            array.close()
        self.arrays.clear()

    @staticmethod
    def pre_load_validation(data_locator):
        location = data_locator.uri_or_path
        if not CxgAdaptor.isvalid(location):
            raise DatasetAccessError(f"cxg matrix is not valid: {location}")

    @staticmethod
    def file_size(data_locator):
        return 0

    @staticmethod
    def open(data_locator, args):
        return CxgAdaptor(data_locator, args)

    def get_about(self):
        return self.about if self.about else super().get_about()

    def get_title(self):
        return self.title if self.title else super().get_title()

    def get_location(self):
        return self.url

    def get_data_locator(self):
        return self.data_locator

    def get_name(self):
        return "cellxgene cxg adaptor version"

    def get_library_versions(self):
        return dict(tiledb=tiledb.__version__)

    def get_path(self, *urls):
        return path_join(self.url, *urls)

    def lsuri(self, uri):
        """
        given a URI, do a tiledb.ls but normalizing for all path weirdness:
        * S3 URIs require trailing slash.  file: doesn't care.
        * results on S3 *have* a trailing slash, Posix does not.

        returns list of (absolute paths, type) *without* trailing slash
        in the path.
        """

        def _cleanpath(p):
            if p[-1] == "/":
                return p[:-1]
            else:
                return p

        if uri[-1] != "/":
            uri += "/"

        result = []
        tiledb.ls(uri, lambda path, type: result.append((_cleanpath(path), type)), ctx=self.tiledb_ctx)
        return result

    @staticmethod
    def isvalid(url):
        """
        Return True if this looks like a valid CXG, False if not.  Just a quick/cheap
        test, not to be fully trusted.
        """
        if not tiledb.object_type(url) == "group":
            return False
        if not tiledb.object_type(path_join(url, "obs")) == "array":
            return False
        if not tiledb.object_type(path_join(url, "var")) == "array":
            return False
        if not tiledb.object_type(path_join(url, "X")) == "array":
            return False
        if not tiledb.object_type(path_join(url, "emb")) == "group":
            return False
        return True

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
        a_type = tiledb.object_type(path_join(self.url, "cxg_group_metadata"))
        if a_type is None:
            # version 0
            cxg_version = "0.0"
            title = None
            about = None
        elif a_type == "array":
            # version >0
            gmd = self.open_array("cxg_group_metadata")
            cxg_version = gmd.meta["cxg_version"]
            if cxg_version == "0.1":
                cxg_properties = json.loads(gmd.meta["cxg_properties"])
                title = cxg_properties.get("title", None)
                about = cxg_properties.get("about", None)

        if cxg_version not in ["0.0", "0.1"]:
            raise DatasetAccessError(f"cxg matrix is not valid: {self.url}")

        self.title = title
        self.about = about
        self.cxg_version = cxg_version

    def open_array(self, name):
        try:
            with self.lock:
                array = self.arrays.get(name)
                if array:
                    return array
                p = self.get_path(name)
                try:
                    array = tiledb.DenseArray(p, mode="r", ctx=self.tiledb_ctx)
                except tiledb.libtiledb.TileDBError as e:
                    raise AttributeError(str(e))
                self.arrays[name] = array
                return array
        except tiledb.libtiledb.TileDBError as e:
            raise AttributeError(str(e))

    def get_embedding_array(self, ename, dims=2):
        array = self.open_array(f"emb/{ename}")
        return array[:, 0:dims]

    def compute_embedding(self, method, filter):
        raise NotImplementedError("CXG does not yet support re-embedding")

    def get_X_array(self, obs_mask=None, var_mask=None):
        obs_items = self._convert_mask(obs_mask)
        var_items = self._convert_mask(var_mask)
        X = self.open_array("X")
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
        except tiledb.libtiledb.TileDBError as e:
            raise AttributeError(str(e))
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

    # function to get the embedding
    # this function to iterate through embeddings.
    def get_embedding_names(self):
        with ServerTiming.time(f"layout.lsuri"):
            pemb = self.get_path("emb")
            embeddings = [os.path.basename(p) for (p, t) in self.lsuri(pemb) if t == "array"]
        return embeddings

    @staticmethod
    def _get_col_type(attr, schema_hints={}):
        type_hint = schema_hints.get(attr.name, {})
        dtype = attr.dtype
        schema = {}
        # type hints take precedence
        if "type" in type_hint:
            schema["type"] = type_hint["type"]
        elif dtype == np.float32:
            schema["type"] = "float32"
        elif dtype == np.int32:
            schema["type"] = "int32"
        elif dtype == np.bool_:
            schema["type"] = "boolean"
        elif dtype == np.str:
            schema["type"] = "string"
        elif dtype == "category":
            schema["type"] = "categorical"
            schema["categories"] = dtype.categories.tolist()
        else:
            raise TypeError(f"Annotations of type {dtype} are unsupported.")

        if schema["type"] == "categorical" and "categories" in schema_hints:
            schema["categories"] = schema_hints["categories"]
        return schema

    def get_schema(self):
        shape = self.get_shape()
        dtype = self.get_X_array_dtype()

        dataframe = {"nObs": shape[0], "nVar": shape[1], "type": dtype.name}

        annotations = {}
        for ax in ("obs", "var"):
            A = self.open_array(ax)
            schema_hints = json.loads(A.meta["cxg_schema"]) if "cxg_schema" in A.meta else {}
            if type(schema_hints) is not dict:
                raise TypeError(f"Array schema was malformed.")

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
                    schema.update(dtype_to_schema(attr.dtype))
                cols.append(schema)

            annotations[ax] = dict(columns=cols)

            if "index" in schema_hints:
                annotations[ax].update({"index": schema_hints["index"]})

        obs_layout = []
        embeddings = self.get_embedding_names()
        for ename in embeddings:
            A = self.open_array(f"emb/{ename}")
            obs_layout.append({"name": ename, "type": A.dtype.name, "dims": [f"{ename}_{d}" for d in range(0, A.ndim)]})

        schema = {"dataframe": dataframe, "annotations": annotations, "layout": {"obs": obs_layout}}
        return schema

    def annotation_to_fbs_matrix(self, axis, fields=None, labels=None):
        with ServerTiming.time(f"annotations.{axis}.query"):
            A = self.open_array(str(axis))
            if axis == Axis.OBS:
                if labels is not None and not labels.empty:
                    df = pd.DataFrame.from_dict(A[:])
                    df = df.join(labels, self.get_obs_names())
                else:
                    df = pd.DataFrame.from_dict(A[:])
            else:
                df = pd.DataFrame.from_dict(A[:])

            if fields is not None and len(fields) > 0:
                df = df[fields]

        with ServerTiming.time(f"annotations.{axis}.encode"):
            fbs = encode_matrix_fbs(df, col_idx=df.columns)

        return fbs

    @staticmethod
    def _convert_mask(boolarray):
        """Convert an index mask to a list of ranges or indices that can be used in a multi_index."""
        if boolarray is None:
            return slice(None)
        assert type(boolarray) == np.ndarray
        assert (boolarray.dtype) == bool

        selector = np.nonzero(boolarray)[0]

        if len(selector) == 0:
            return slice(None)

        result = []
        current = slice(selector[0], selector[0])
        for sel in selector[1:]:
            if sel == current.stop + 1:
                current = slice(current.start, sel)
            else:
                result.append(current if current.start != current.stop else current.start)
                current = slice(sel, sel)

        if len(result) == 0 or result[-1] != current:
            result.append(current if current.start != current.stop else current.start)

        return result
