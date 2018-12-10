import click

from .launch import launch
from .prepare import prepare


@click.group(name="cellxgene", context_settings=dict(max_content_width=85))
@click.version_option(version="0.3.0", prog_name="cellxgene", message="[%(prog)s] Version %(version)s")
def cli():
    pass


cli.add_command(launch)
cli.add_command(prepare)
