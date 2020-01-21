import contextlib
import errno
import socket
import copy
from server.data_common.fbs.matrix import decode_matrix_fbs
from server.common.errors import DisabledFeatureError


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


def get_schema(data, annotations=None, session=None):
    """helper function to gather the schema from the data source and annotations"""
    schema = data.get_schema()
    # add label obs annotations as needed
    if annotations is not None:
        labels = annotations.read_labels(session)
        if labels is not None and not labels.empty:
            schema = copy.deepcopy(schema)
            for col in labels.columns:
                col_schema = {
                    "name": col,
                    "writable": True,
                }
                # FIXME: data._get_col_type needs to be refactored.
                col_schema.update(data._get_col_type(labels[col]))
                schema["annotations"]["obs"]["columns"].append(col_schema)

    return schema


def annotation_put_fbs(fbs, data, annotations=None, session=None):
    """helper function to write annotations from fbs"""
    if annotations is None:
        raise DisabledFeatureError("Writable annotations are not enabled")

    new_label_df = decode_matrix_fbs(fbs)
    if not new_label_df.empty:
        # FIXME - refactor original_obs_index
        new_label_df.index = data.original_obs_index

    data.validate_label_data(new_label_df)
    annotations.write_labels(new_label_df, session)
