from functools import wraps

from flask import json
from numpy import float32, integer

from server.common.errors import DataAdaptorError
from enum import Enum


class Float32JSONEncoder(json.JSONEncoder):
    def __init__(self, *args, **kwargs):
        """
        NaN/Infinities are illegal in standard JSON.  Python extends JSON with
        non-standard symbols that most JavaScript JSON parsers do not understand.
        The `allow_nan` parameter will force Python simplejson to throw an ValueError
        if it runs into non-finite floating point values which are unsupported by
        standard JSON.
        """
        kwargs["allow_nan"] = False
        super().__init__(*args, **kwargs)

    def default(self, obj):
        if isinstance(obj, float32):
            return float(obj)
        elif isinstance(obj, integer):
            return int(obj)
        return json.JSONEncoder.default(self, obj)


def custom_format_warning(msg, *args, **kwargs):
    return f"[cellxgene] Warning: {msg} \n"


def jsonify_numpy(data):
    return json.dumps(data, cls=Float32JSONEncoder, allow_nan=False)


class MatrixDataType(Enum):
    H5AD = "h5ad"
    TILEDB = "tiledb"
    UNKNOWN = "unknown"


class MatrixDataLoader(object):

    def __init__(self, location):
        self.location = location
        self.etype = self.matrix_data_type()
        self.matrix_type = None
        if self.etype == MatrixDataType.H5AD:
            from server.data_scanpy.scanpy_adaptor import ScanpyAdaptor
            self.matrix_type = ScanpyAdaptor
        elif self.etype == MatrixDataType.TILEDB:
            from server.data_tiledb.tiledb_adaptor import TileDbAdaptor
            self.matrix_type = TileDbAdaptor

    def matrix_data_type(self):
        if self.location.endswith(".h5ad"):
            return MatrixDataType.H5AD
        elif ".cxg" in self.location:
            return MatrixDataType.TILEDB
        elif self.location.startswith("s3://"):
            return MatrixDataType.TILEDB
        else:
            return MatrixDataType.UNKNOWN

    def pre_checks(self):
        if self.etype == MatrixDataType.UNKNOWN:
            raise RuntimeError(f"{self.location} does not have a recognized type: .h5ad or .cxg")
        self.matrix_type.pre_checks(self.location)

    def file_size(self):
        return self.matrix_type.file_size(self.location)

    def open(self, app_config):
        try:
            # create and return an object to the matrix_type
            return self.matrix_type.open(self.location, app_config)
        except Exception as e:
            raise RuntimeError(str(e))
