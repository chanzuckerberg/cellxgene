import contextlib
import errno
import socket

def find_available_port(host, port):
    """
    Helper method to find open port on host. Tries 50 ports incremented from the specified port
    """
    new_port = port
    found_port = False
    num_ports_to_try = 50
    with contextlib.closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as s:
        for port_to_try in range(port, port + num_ports_to_try):
            while True:
                try:
                    s.bind((host, port_to_try))
                    found_port = True
                    new_port = port_to_try
                    break
                except socket.error as e:
                    # Reraise if error is not "address in use"
                    if e.errno != errno.EADDRINUSE:
                        raise e
                    break
            if found_port:
                break
    if not found_port:
        print(f"[cellxgene] No port in range {new_port} - {new_port + num_ports_to_try - 1} was"
                   f" available to run cellxgene server on.")
        print(f"[cellxgene] Exiting")
        raise OSError("No available ports")
    return new_port
