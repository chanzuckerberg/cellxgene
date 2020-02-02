import os
import json
from server.data_common.driver import CXGDriver
import tiledb
import numpy as np
import pandas as pd
from server.data_common.fbs.matrix import encode_matrix_fbs
from server.common.utils import path_join
from server_timing import Timing as ServerTiming
from server.common.constants import Axis


class TileDbEngine(CXGDriver):
    def __init__(self, location=None, config=None):
        super().__init__(config)
        self.url = location
        self.update(location, config)
        self.tiledb_ctx = tiledb.Ctx({
            'sm.tile_cache_size': 8 * 1024 * 1024 * 1024,
            'sm.num_reader_threads': 32,
        })

    def update(self, location=None, config=None):
        self.config = config
        if location:
            self.url = location
            if self.url[-1] != '/':
                self.url += '/'

            self._validate_and_initialize()

    @staticmethod
    def pre_checks(location):
        if not TileDbEngine.isvalid(location):
            raise RuntimeError(f"tiledb matrix is not valid: {location}")

    @staticmethod
    def file_size(location):
        return 0

    @staticmethod
    def open(location, args):
        return TileDbEngine(location, args)

    def get_location(self):
        return self.url

    def get_name(self):
        return "cellxgene TileDb engine version"

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
            if p[-1] == '/':
                return p[:-1]
            else:
                return p

        if uri[-1] != '/':
            uri += '/'

        result = []
        tiledb.ls(uri,
                  lambda path, type: result.append((_cleanpath(path), type)),
                  ctx=self.tiledb_ctx)
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
        if not self.isvalid(self.url):
            raise RuntimeError(f"invalid tiledb dataset {self.url}")

    def get_array(self, name, items=None):
        p = self.get_path(name)
        if items is None:
            items = slice(None)
        with tiledb.DenseArray(p, mode="r", ctx=self.tiledb_ctx) as A:
            data = A[items]
        return data

    def get_embedding_array(self, ename, items=None):
        return self.get_array(f"emb/{ename}", items)

    def get_X_array(self, obs_items=None, var_items=None):
        p = self.get_path("X")
        if obs_items is None:
            obs_items = slice(None)
        if var_items is None:
            var_items = slice(None)
        # FIXME, the is a problem retrieving index sets from tiledb
        with tiledb.DenseArray(p, mode="r", ctx=self.tiledb_ctx) as X:
            tdata = X[:, :]
            data = tdata[obs_items, var_items]
        return data

    def get_X_array_shape(self):
        p = self.get_path("X")
        with tiledb.DenseArray(p, mode="r", ctx=self.tiledb_ctx) as A:
            return A.shape

    def get_X_array_dtype(self):
        p = self.get_path("X")
        with tiledb.DenseArray(p, mode="r", ctx=self.tiledb_ctx) as A:
            return A.dtype

    def query_var_array(self, term_name):
        p = self.get_path("var")
        try:
            with tiledb.DenseArray(p, mode='r', ctx=self.tiledb_ctx) as A:
                anno_data = A.query(attrs=[term_name])[:][term_name]
        except tiledb.libtiledb.TileDBError as e:
            raise AttributeError(str(e))
        return anno_data

    def query_obs_array(self, term_name):
        p = self.get_path("obs")
        try:
            with tiledb.DenseArray(p, mode='r', ctx=self.tiledb_ctx) as A:
                anno_data = A.query(attrs=[term_name])[:][term_name]
        except tiledb.libtiledb.TileDBError as e:
            raise AttributeError(str(e))
        return anno_data

    def get_obs_names(self):
        # get the index from the meta data
        p = self.get_path("obs")
        with tiledb.DenseArray(p, mode='r', ctx=self.tiledb_ctx) as A:
            meta = json.loads(A.meta["cxg_schema"])
            index_name = meta["index"]
            return index_name

    def get_obs_index(self):
        p = self.get_path("obs")
        with tiledb.DenseArray(p, mode='r', ctx=self.tiledb_ctx) as A:
            meta = json.loads(A.meta["cxg_schema"])
            index_name = meta["index"]
            data = A.query(attrs=[index_name])[:][index_name]
            return data

    def get_obs_columns(self):
        p = self.get_path("obs")
        with tiledb.DenseArray(p, mode='r', ctx=self.tiledb_ctx) as A:
            schema = A.schema
            col_names = [attr.name for attr in schema]
            return pd.Index(col_names)

    def get_obs_shape(self):
        p = self.get_path("obs")
        with tiledb.DenseArray(p, mode='r', ctx=self.tiledb_ctx) as A:
            return A.shape

    # function to get the embedding
    # this function to iterate through embeddings.
    def get_embedding_names(self):
        with ServerTiming.time(f'layout.lsuri'):
            pemb = self.get_path("emb")
            embeddings = [
                os.path.basename(p) for (p, t) in self.lsuri(pemb) if t == 'array'
            ]
        return embeddings

    @staticmethod
    def _get_col_type(attr, schema_hints={}):
        type_hint = schema_hints.get(attr.name, {})
        dtype = attr.dtype
        schema = {}
        # type hints take precedence
        if 'type' in type_hint:
            schema['type'] = type_hint['type']
        elif dtype == np.float32:
            schema['type'] = 'float32'
        elif dtype == np.int32:
            schema['type'] = 'int32'
        elif dtype == np.bool_:
            schema['type'] = 'boolean'
        elif dtype == np.str:
            schema['type'] = 'string'
        elif dtype == "category":
            schema["type"] = "categorical"
            schema["categories"] = dtype.categories.tolist()
        else:
            raise TypeError(
                f"Annotations of type {dtype} are unsupported."
            )

        if schema['type'] == 'categorical' and 'categories' in schema_hints:
            schema['categories'] = schema_hints['categories']
        return schema

    def get_schema(self):
        shape = self.get_X_array_shape()
        dtype = self.get_X_array_dtype()

        dataframe = {
            'nObs': shape[0],
            'nVar': shape[1],
            'type': dtype.name
        }

        annotations = {}
        for ax in ('obs', 'var'):
            p = self.get_path(ax)
            with tiledb.Array(p, ctx=self.tiledb_ctx) as A:
                schema_hints = json.loads(A.meta['cxg_schema']) if 'cxg_schema' in A.meta else {}
                if type(schema_hints) is not dict:
                    raise TypeError(f'Array schema was malformed.')

                annotations[ax] = {
                    'columns': [
                        {
                            'name': attr.name,
                            'writable': False,
                            **self._get_col_type(attr, schema_hints)
                        } for attr in A.schema
                    ]
                }
                if 'index' in schema_hints:
                    annotations[ax].update({'index': schema_hints['index']})

        obs_layout = []
        embeddings = self.get_embedding_names()
        for ename in embeddings:
            with tiledb.DenseArray(self.get_path(f"emb/{ename}"), ctx=self.tiledb_ctx) as A:
                obs_layout.append({
                    'name': ename,
                    'type': A.dtype.name,
                    'dims': [f'{ename}_{d}' for d in range(0, A.ndim)]
                })

        schema = {
            'dataframe': dataframe,
            'annotations': annotations,
            'layout': {'obs': obs_layout}
        }
        return schema

    def annotation_to_fbs_matrix(self, axis, fields=None, labels=None):
        with ServerTiming.time(f'annotations.{axis}.query'):
            p = self.get_path(str(axis))
            with tiledb.DenseArray(p, mode='r', ctx=self.tiledb_ctx) as A:
                if fields is not None and len(fields) > 0:
                    try:
                        df = pd.DataFrame(A.query(attrs=fields)[:])
                    except tiledb.libtiledb.TileDBError:
                        raise KeyError("bad field {fields}")

                else:
                    df = pd.DataFrame.from_dict(A[:])

                if axis == Axis.OBS:
                    if labels is not None and not labels.empty:
                        obs_names = self.get_obs_names()
                        df = df.join(labels, obs_names)

        with ServerTiming.time(f'annotations.{axis}.encode'):
            fbs = encode_matrix_fbs(df, col_idx=df.columns)

        return fbs
