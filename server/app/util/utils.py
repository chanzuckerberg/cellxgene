import json
from argparse import ArgumentTypeError

from numpy import float32, integer

from server.app.util.errors import MimeTypeError


class Float32JSONEncoder(json.JSONEncoder):
    def __init__(self, *args, **kwargs):
        """
        NaN/Infinities are illegal in standard JSON.  Python extends JSON with
        non-standard symbols that most JavaScript JSON parsers do not understand.
        The `allow_nan` parameter will force Python simplejson to throw an ValueError
        if it runs into non-finite floating point values which are unsupported by
        standard JSON.
        """
        kwargs['allow_nan'] = False
        super().__init__(*args, **kwargs)

    def default(self, obj):
        if isinstance(obj, float32):
            return float(obj)
        elif isinstance(obj, integer):
            return int(obj)
        return json.JSONEncoder.default(self, obj)


def custom_format_warning(msg, *args, **kwargs):
    return f"[cellxgene] Warning: {msg} \n"


def get_mime_type(default="application/json", acceptable_types=["application/json", "text/csv"], query_param=None,
                  header=None):
    mime_type = default
    if query_param:
        if query_param in acceptable_types:
            mime_type = query_param
        else:
            raise MimeTypeError(f"Unsupported mime type {query_param} specified in query parameter 'accept-type'")
    elif len(header):
        mime_type = header.best_match(acceptable_types)
        if not mime_type:
            raise MimeTypeError(f"Unsupported mime type(s) {header} in HTTP Accept header")
    return mime_type


def whole_number(value):
    try:
        value = int(value)
    except ValueError as e:
        raise ArgumentTypeError(f"{value} is not type int") from e
    if value < 0:
        raise ArgumentTypeError(f"{value} is not >= 0")
    return value
