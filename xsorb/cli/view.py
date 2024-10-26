'''
CLI parser for command: view
'''

import argparse

from xsorb.cli.command import CLICommandBase, nonnegative_int


class CLICommand(CLICommandBase):
    """open the selected calculation with ASE GUI

Example:
 $ xsorb view relax 12
    """

    @staticmethod
    def add_arguments(parser : argparse.ArgumentParser):
        parser.add_argument('calc_type',
                        type=str,
                        choices=['screening', 'relax', 'ml_opt', 'initial'],
                        help='''type of calculation to plot''')
        parser.add_argument('calc_id',
                        nargs='?',
                        type=nonnegative_int,
                        help='''ID of the calculation to plot''')

    @staticmethod
    def run(args : argparse.Namespace):
        from xsorb.visualize.images import view_config
        view_config(args.calc_type, args.calc_id)


    @staticmethod
    def bind_function(parser: argparse.ArgumentParser):
        parser.set_defaults(func=CLICommand.run)
