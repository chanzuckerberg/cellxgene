import click


def find_available_port(host, port):
    """
    Helper method to find open port on host. Tries 50 ports incremented from the specified port
    """
    import errno
    import socket
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    new_port = port
    found_port = False
    num_ports_to_try = 50
    for port_to_try in range(port, port + num_ports_to_try):
        new_port = port + i
        while True:
            try:
                s.bind((host, new_port))
                found_port = True
                break
            except socket.error as e:
                # Reraise if error is not "address in use"
                if e.errno != errno.EADDRINUSE:
                    raise e
                break
        if found_port:
            break
    if not found_port:
        click.echo(f"[cellxgene] No port in range {new_port} - {new_port + num_ports_to_try - 1} was"
                   f" available to run cellxgene server on.")
        click.echo(f"[cellxgene] Exiting")
        raise OSError("No available ports")
    s.close()
    return new_port
