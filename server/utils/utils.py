import contextlib
import errno
import socket


def find_available_port(host, port=5005):
    """
    Helper method to find open port on host. Tries 5000 ports incremented from the specified port
    """
    # Takes approx 2 seconds to do a scan of 5000 ports on my laptop
    num_ports_to_try = 5000
    for port_to_try in range(port, port + num_ports_to_try):
        with contextlib.closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as s:
            try:
                s.bind((host, port_to_try))
                return port_to_try
            except socket.error:
                pass
    raise OSError(errno.EADDRINUSE, f"No port in range {port} - {port + num_ports_to_try - 1} available.")
