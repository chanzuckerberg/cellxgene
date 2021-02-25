import click

from local_server.converters.schema import remix, validate


@click.group(
    name="schema",
    subcommand_metavar="COMMAND <args>",
    short_help="Apply and validate the cellxgene data integration schema to an h5ad file.",
    context_settings=dict(max_content_width=85, help_option_names=["-h", "--help"]),
)
def schema_cli():
    try:
        import scanpy  # noqa: F401
    except ImportError:
        raise click.ClickException(
            "[cellxgene] cellxgene schema requires scanpy"
        )


@click.command(
    name="apply",
    short_help="(experimental) Apply the cellxgene data integration schema to an h5ad.",
    help="(experimental) Using a yaml file that describes schema values to insert or convert and in input "
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
    name="validate",
    short_help="(experimental) Check that an h5ad follows the cellxgene data integration schema.",
)
@click.argument(
    "h5ad",
    nargs=1,
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "--shallow",
    help="When true, just check that the correct version information is present.",
    default=False,
    show_default=True,
    is_flag=True,
)
def schema_validate(h5ad, shallow):
    validate.validate(h5ad, shallow)


schema_cli.add_command(schema_apply)
schema_cli.add_command(schema_validate)
