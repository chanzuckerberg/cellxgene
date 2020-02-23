import contextlib
import errno
import socket
from urllib.parse import urlsplit, urljoin
import os
from flask import json
import numpy as np
import pandas as pd
import warnings


def find_available_port(host, port=5005):
    """
    Helper method to find open port on host. Tries 5000 ports incremented from the specified port
    """
    # Takes approx 2 seconds to do a scan of 5000 ports on my laptop
    num_ports_to_try = 5000
    for port_to_try in range(port, port + num_ports_to_try):
        if is_port_available(host, port_to_try):
            return port_to_try
    raise socket.error(errno.EADDRINUSE, f"No port in range {port} - {port + num_ports_to_try - 1} available.")


def is_port_available(host, port):
    is_available = False
    with contextlib.closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as s:
        try:
            s.bind((host, port))
            is_available = True
        except socket.error:
            pass
    return is_available


def sort_options(command):
    """
    Helper for the click options - will sort options in a command, and can
    be used as a decorator.
    """
    command.params.sort(key=lambda p: p.name)
    return command


def path_join(base, *urls):
    """
    this is like urllib.parse.urljoin, except it works around the scheme-specific
    cleverness in the aforementioned code, ignores anything in the url except the path,
    and accepts more than one url.
    """
    if not base.endswith("/"):
        base += "/"
    btpl = urlsplit(base)
    path = btpl.path
    for url in urls:
        utpl = urlsplit(url)
        if btpl.scheme == "":
            path = os.path.join(path, utpl.path)
            path = os.path.normpath(path)
        else:
            path = urljoin(path, utpl.path)
    return btpl._replace(path=path).geturl()


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
        if isinstance(obj, np.float32):
            return float(obj)
        elif isinstance(obj, np.integer):
            return int(obj)
        return json.JSONEncoder.default(self, obj)


def custom_format_warning(msg, *args, **kwargs):
    return f"[cellxgene] Warning: {msg} \n"


def jsonify_numpy(data):
    return json.dumps(data, cls=Float32JSONEncoder, allow_nan=False)


def dtype_to_schema(dtype):
    schema = {}
    if dtype == np.float32:
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
    return schema


def can_cast_to_float32(array):
    if array.dtype.kind == "f":
        if not np.can_cast(array.dtype, np.float32):
            warnings.warn(f"Annotation {array.name} will be converted to 32 bit float and may lose precision.")
        return True
    return False


def can_cast_to_int32(array):
    if array.dtype.kind in ["i", "u"]:
        if np.can_cast(array.dtype, np.int32):
            return True
        ii32 = np.iinfo(np.int32)
        if array.min() >= ii32.min and array.max() <= ii32.max:
            return True
    return False


def series_to_schema(array):
    assert type(array) == pd.Series
    try:
        return dtype_to_schema(array.dtype)
    except TypeError:
        dtype = array.dtype
        data_kind = dtype.kind
        schema = {}
        if can_cast_to_float32(array):
            schema["type"] = "float32"
        elif can_cast_to_int32(array):
            schema["type"] = "int32"
        elif data_kind == "O" and dtype == "object":
            schema["type"] = "string"
        else:
            raise TypeError(f"Annotations of type {dtype} are unsupported.")
        return schema
