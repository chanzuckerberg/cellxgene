import click


@click.group(
    name="schema",
    subcommand_metavar="COMMAND <args>",
    short_help="Apply and validate teh cellxgene data integration schema to an h5ad file.",
    context_settings=dict(max_content_width=85, help_option_names=["-h", "--help"]),
)
def schema_cli():
    pass


@click.command(
    name="apply",
    short_help="Apply the cellxgene data integration schema to an h5ad.",
    help="Using a yaml file that describes schema values to insert or convert and in input "
    "h5ad file, apply the schema changes and create a new, conforming h5ad.",
)
@click.argument(
    "--source-h5ad",
    nargs=1,
    type=click.Path(exists=True, dir_okay=False),
)
@click.argument(
    "--remix-config",
    nargs=1,
    type=click.Path(exists=True, dir_okay=False),
)
@click.argument("--output_filename", nargs=1)
def schema_apply(source_h5ad, remix_config, output_filename):
    pass


@click.command(
    name="check",
    short_help="Check that an h5ad follows the cellxgene data integration schema.",
)
def schema_check():
    pass


schema_cli.add_command(schema_apply)
schema_cli.add_command(schema_check)
