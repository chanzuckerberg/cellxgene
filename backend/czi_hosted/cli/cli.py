import click

from .convert_to_cxg import convert_to_cxg
from .launch import launch
from .prepare import prepare
from .upgrade import log_upgrade_check
from .schema import schema_cli
from .. import __version__


@click.group(
    name="cellxgene",
    subcommand_metavar="COMMAND <args>",
    options_metavar="<options>",
    context_settings=dict(max_content_width=85, help_option_names=["-h", "--help"]),
)
@click.help_option("--help", "-h", help="Show this message and exit.")
@click.version_option(
    version=__version__,
    prog_name="cellxgene",
    message="[%(prog)s] Version %(version)s",
    help="Show the software version and exit.",
)
@click.option(
    "--upgrade-check/--no-upgrade-check", default=True, show_default=True, help="Check for release upgrades on start.",
)
def cli(upgrade_check):
    if upgrade_check:
        log_upgrade_check()


cli.add_command(launch)
cli.add_command(prepare)
cli.add_command(convert_to_cxg)
cli.add_command(schema_cli)
