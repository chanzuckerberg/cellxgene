import contextlib
import errno
import socket
from urllib.parse import urlsplit, urljoin
import os
from flask import json
from numpy import float32, integer


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
    btpl = urlsplit(base)
    path = btpl.path
    for url in urls:
        utpl = urlsplit(url)
        if btpl.scheme == "":
            path = os.path.join(path, utpl.path)
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
        if isinstance(obj, float32):
            return float(obj)
        elif isinstance(obj, integer):
            return int(obj)
        return json.JSONEncoder.default(self, obj)


def custom_format_warning(msg, *args, **kwargs):
    return f"[cellxgene] Warning: {msg} \n"


def jsonify_numpy(data):
    return json.dumps(data, cls=Float32JSONEncoder, allow_nan=False)
