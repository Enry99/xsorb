'''
CLI parser for command: status
'''

import argparse

from xsorb.cli.command import CLICommandBase


class CLICommand(CLICommandBase):
    """get information about the status of calculations

Example:
 $ xsorb status
    """

    @staticmethod
    def add_arguments(parser : argparse.ArgumentParser):
        parser.add_argument('calc_type',
                        type=str,
                        nargs='?',
                        default='all',
                        choices=['screening', 'relax', 'mlopt', 'all'],
                        help='''type of calculation to plot''')

    @staticmethod
    def run(args : argparse.Namespace):
        from xsorb.io.status import print_status # pylint: disable=import-outside-toplevel
        print_status(args.calc_type)


    @staticmethod
    def bind_function(parser: argparse.ArgumentParser):
        parser.set_defaults(func=CLICommand.run)
