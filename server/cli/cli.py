import click

from .launch import launch
from .prepare import prepare


@click.group(name="cellxgene",
             subcommand_metavar="COMMAND <args>",
             options_metavar="<options>",
             context_settings=dict(max_content_width=85,
                                   help_option_names=['-h', '--help']))
@click.help_option("--help", "-h", help="Show this message and exit.")
@click.version_option(
    version="0.13.0",
    prog_name="cellxgene",
    message="[%(prog)s] Version %(version)s",
    help="Show the software version and exit.")
def cli():
    pass


cli.add_command(launch)
cli.add_command(prepare)
