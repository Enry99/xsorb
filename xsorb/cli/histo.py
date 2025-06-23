'''
CLI parser for command: histo
'''

import argparse

from xsorb.cli.command import CLICommandBase


class CLICommand(CLICommandBase):
    """save image with adsorption energy histogram for a given calculation type
Example:
 $ xsorb histo screening
    """

    @staticmethod
    def add_arguments(parser : argparse.ArgumentParser):
        parser.add_argument('calc_type',
                        type=str,
                        choices=['screening', 'relax', 'mlopt'],
                        help='''type of calculation''')


    @staticmethod
    def run(args : argparse.Namespace):
        from xsorb.visualize.images import plot_histo # pylint: disable=import-outside-toplevel
        plot_histo(args.calc_type)


    @staticmethod
    def bind_function(parser: argparse.ArgumentParser):
        parser.set_defaults(func=CLICommand.run)
