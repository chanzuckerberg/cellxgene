from functools import wraps

from flask import json
from numpy import float32, integer

from server.app.util.errors import DriverError


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


def jsonify_scanpy(data):
    return json.dumps(data, cls=Float32JSONEncoder, allow_nan=False)


def requires_data(func):
    @wraps(func)
    def wrapped_function(self, *args, **kwargs):
        if self.data is None:
            raise DriverError(f"error data must be loaded before you call {func.__name__}")
        return func(self, *args, **kwargs)
    return wrapped_function
