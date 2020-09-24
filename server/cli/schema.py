import click

from server.converters.schema import remix


@click.group(
    name="schema",
    subcommand_metavar="COMMAND <args>",
    short_help="Apply and validate the cellxgene data integration schema to an h5ad file.",
    context_settings=dict(max_content_width=85, help_option_names=["-h", "--help"]),
)
def schema_cli():
    try:
        import scanpy
    except ImportError:
        raise click.ClickException(
            "[cellxgene] cellxgene schema requires scanpy"
        )
    scanpy.__version__


@click.command(
    name="apply",
    short_help="Apply the cellxgene data integration schema to an h5ad.",
    help="Using a yaml file that describes schema values to insert or convert and in input "
    "h5ad file, apply the schema changes and create a new, conforming h5ad.",
)
@click.option(
    "--source-h5ad",
    help="Input h5ad file.",
    nargs=1,
    required=True,
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "--remix-config",
    help="Config yaml with information on how to apply the schema.",
    nargs=1,
    required=True,
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "--output-filename",
    help="Filename for the new, schema-conforming h5ad file.",
    required=True,
    nargs=1
)
def schema_apply(source_h5ad, remix_config, output_filename):
    remix.apply_schema(source_h5ad, remix_config, output_filename)


@click.command(
    name="check",
    short_help="Check that an h5ad follows the cellxgene data integration schema.",
)
def schema_check():
    pass


schema_cli.add_command(schema_apply)
schema_cli.add_command(schema_check)
