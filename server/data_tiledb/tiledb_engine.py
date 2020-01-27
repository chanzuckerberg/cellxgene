import os
import json
from server.data_common.driver import CXGDriver
from urllib.parse import urlsplit, urljoin
import tiledb
from server_timing import Timing as ServerTiming
import numpy as np
import pandas as pd
from server.data_common.fbs.matrix import encode_matrix_fbs


class TileDbEngine(CXGDriver):
    def __init__(self, location=None, args={}):
        self.url = location
        self.update(location, args)
        self.config = self._get_default_config()
        self.config.update(args)
        self.tiledb_ctx = tiledb.Ctx({
            'sm.tile_cache_size': 8 * 1024 * 1024 * 1024,
            'sm.num_reader_threads': 32,
        })

    def update(self, location=None, args={}):
        if location:
            self.url = location
            if self.url[-1] != '/':
                self.url += '/'

            self._validate_and_initialize()

    @staticmethod
    def pre_check(location):
        pass

    @staticmethod
    def file_size(location):
        return 0

    @staticmethod
    def open(location, args):
        return TileDbEngine(location, args)

    @staticmethod
    def _get_default_config():
        return {
            "layout": [],
            "max_category_items": 100,
            "obs_names": None,
            "var_names": None,
            "diffexp_lfc_cutoff": 0.01,
            "backed": False,
            "disable_diffexp": False,
            "diffexp_may_be_slow": False,
        }

    @staticmethod
    def join(base, *urls):
        """
        this is like urllib.parse.urljoin, except it works around the scheme-specific
        cleverness in the aforementioned code, ignores anything in the url except the path,
        and accepts more than one url.
        """
        btpl = urlsplit(base)
        path = btpl.path
        for url in urls:
            utpl = urlsplit(url)
            path = urljoin(path, utpl.path)
        return btpl._replace(path=path).geturl()

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

    def isvalid(self):
        """
        Return True if this looks like a valid CXG, False if not.  Just a quick/cheap
        test, not to be fully trusted.
        """
        if not tiledb.object_type(self.url) == "group":
            return False
        if not tiledb.object_type(self.obs) == "array":
            return False
        if not tiledb.object_type(self.var) == "array":
            return False
        if not tiledb.object_type(self.X) == "array":
            return False
        if not tiledb.object_type(self.emb) == "group":
            return False
        return True

    def _validate_and_initialize(self):
        pass

    def get_config_parameters(self, uid=None):
        # TODO
        params = {
            "max-category-items": self.config["max_category_items"],
            "disable-diffexp": self.config["disable_diffexp"],
            "diffexp-may-be-slow": self.config["diffexp_may_be_slow"],
        }
        return params

    @property
    def obs(self):
        return self.join(self.url, 'obs')

    @property
    def var(self):
        return self.join(self.url, 'var')

    @property
    def X(self):
        return self.join(self.url, 'X')

    @property
    def emb(self):
        return self.join(self.url, 'emb')

    def embedding(self, ename):
        return self.join(self.url, 'emb/', ename)

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
        else:
            raise TypeError(
                f"Annotations of type {dtype} are unsupported."
            )

        if schema['type'] == 'categorical' and 'categories' in schema_hints:
            schema['categories'] = schema_hints['categories']
        return schema

    def get_schema(self, uid=None, collection=None):
        with tiledb.Array(self.X, ctx=self.tiledb_ctx) as X:
            dataframe = {
                'nObs': X.shape[0],
                'nVar': X.shape[1],
                'type': str(X.dtype)
            }

        annotations = {}
        for ax in ('obs', 'var'):
            with tiledb.Array(getattr(self, ax), ctx=self.tiledb_ctx) as A:
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
        embeddings = [
            os.path.basename(p) for (p, t) in self.lsuri(self.emb) if t == 'array'
        ]
        for ename in embeddings:
            with tiledb.DenseArray(self.embedding(ename), ctx=self.tiledb_ctx) as A:
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

    def layout_to_fbs_matrix(self):
        """
        return all embeddings as a flatbuffer, using the cellxgene matrix fbs encoding.

        * returns only first two dimensions, with name {ename}_0 and {ename}_1,
          where {ename} is the embedding name.
        * client assumes each will be individually centered & scaled (isotropically)
          to a [0, 1] range.

        """
        with ServerTiming.time(f'layout.lsuri'):
            embeddings = [
                os.path.basename(p) for (p, t) in self.lsuri(self.emb) if t == 'array'
            ]

        layout_data = []
        with ServerTiming.time(f'layout.query'):
            for ename in embeddings:
                with tiledb.DenseArray(self.embedding(ename), mode='r', ctx=self.tiledb_ctx) as E:
                    embedding = E[:, 0:2]

                    # scale isotropically
                    min = embedding.min(axis=0)
                    max = embedding.max(axis=0)
                    scale = np.amax(max - min)
                    normalized_layout = (embedding - min) / scale

                    # translate to center on both axis
                    translate = 0.5 - ((max - min) / scale / 2)
                    normalized_layout = normalized_layout + translate

                    normalized_layout = normalized_layout.astype(dtype=np.float32)
                    layout_data.append(pd.DataFrame(normalized_layout, columns=[f"{ename}_0", f"{ename}_1"]))

        with ServerTiming.time(f'layout.encode'):
            if layout_data:
                df = pd.concat(layout_data, axis=1, copy=False)
            else:
                df = pd.DataFrame()
            fbs = encode_matrix_fbs(df, col_idx=df.columns, row_idx=None)

        return fbs

    def annotation_to_fbs_matrix(self, axis, fields=None, labels=None):
        with ServerTiming.time(f'annotations.{axis}.query'):
            with tiledb.DenseArray(getattr(self, str(axis)), mode='r', ctx=self.tiledb_ctx) as A:
                # FIXME handle labels
                if not fields:
                    data = A[:]
                else:
                    data = A.query(attrs=fields)[:]

                df = pd.DataFrame.from_dict(data)

        with ServerTiming.time(f'annotations.{axis}.encode'):
            fbs = encode_matrix_fbs(df, col_idx=df.columns)

        return fbs

    def data_frame_to_fbs_matrix(self, filter, axis):
        if axis != 'var':
            # not yet implemented in the FBS encoder.
            raise ValueError('Obs dimension slicing not yet supported')

        with ServerTiming.time(f'data.{axis}.query'):
            with tiledb.DenseArray(self.X, mode='r', ctx=self.tiledb_ctx) as X:
                shape = X.shape
                (obs_selector, var_selector) = self.filter_parse(filter, shape)

                # multi_index only accepts lists or slices
                oidx = obs_selector.tolist() if type(obs_selector) is np.ndarray else obs_selector
                vidx = var_selector.tolist() if type(var_selector) is np.ndarray else var_selector
                data = X.multi_index[oidx, vidx]['']

        with ServerTiming.time(f'data.{axis}.encode'):
            fbs = encode_matrix_fbs(data, col_idx=var_selector, row_idx=None)

        return fbs

    def diffexp_topN(self, *args):
        # TODO
        raise RuntimeError("Not implemented yet")
