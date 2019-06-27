import shlex

import click

from server.cli.launch import common_args
from server.gui.utils import OptionsError


@click.command()
@common_args
def cli(**kwargs):
    pass


def parse_opt_string(opts):
    context = click.Context(cli)
    parser = click.OptionParser(context)
    for command in context.command.params:
        command.add_to_parser(parser, context)

    try:
        opts, args, param_order = parser.parse_args(shlex.split(opts))
        for param in cli.params:
            value, args = param.handle_parse_result(context, opts, args)
    except click.ClickException as ce:
        raise OptionsError(ce.message) from ce
    return dict(context.params)
