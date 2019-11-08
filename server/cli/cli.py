import click

from .launch import launch
from .prepare import prepare


@click.group(name="cellxgene",
             context_settings=dict(max_content_width=85,
                                   help_option_names=['-h', '--help']))
@click.help_option("--help", "-h", help="Show this message and exit.")
@click.version_option(
    version="0.12.0",
    prog_name="cellxgene",
    message="[%(prog)s] Version %(version)s",
    help="Show the software version and exit.")
def cli():
    pass


cli.add_command(launch)
cli.add_command(prepare)
