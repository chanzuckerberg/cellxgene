import contextlib
import errno
import socket


def find_available_port(host, port):
    """
    Helper method to find open port on host. Tries 50 ports incremented from the specified port
    """
    num_ports_to_try = 50
    for port_to_try in range(port, port + num_ports_to_try):
        with contextlib.closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as s:
            while True:
                try:
                    s.bind((host, port_to_try))
                    return port_to_try
                except socket.error as e:
                    # Reraise if error is not "address in use"
                    if e.errno != errno.EADDRINUSE:
                        raise e
                    break
    print(f"[cellxgene] No port in range {port} - {port + num_ports_to_try - 1} was"
          f" available to run cellxgene server on.")
    print(f"[cellxgene] Exiting")
    raise OSError("No available ports")
