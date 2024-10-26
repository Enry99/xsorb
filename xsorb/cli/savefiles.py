'''
CLI parser for command: savefiles
'''

import argparse

from xsorb.cli.command import CLICommandBase


class CLICommand(CLICommandBase):
    """save files in specific format

Example:
 $ xsorb savefiles screening cif
    """

    @staticmethod
    def add_arguments(parser : argparse.ArgumentParser):
        parser.add_argument('calc_type',
                        type=str,
                        choices=['screening', 'relax', 'ml_opt', 'initial'],
                        help='''type of calculation to render''')
        parser.add_argument('format',
                        type=str,
                        help='''format to save the files''')

    @staticmethod
    def run(args : argparse.Namespace):
        from xsorb.io.inputs import saveas
        saveas(args.calc_type, saveas_format=args.format)


    @staticmethod
    def bind_function(parser: argparse.ArgumentParser):
        parser.set_defaults(func=CLICommand.run)
