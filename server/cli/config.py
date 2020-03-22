import click
from server.common.default_config import default_config


@click.command(short_help="Print the default cellxgene server configuration")
def config():
    print(default_config)
