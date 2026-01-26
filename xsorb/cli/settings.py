'''
CLI parser for command: status
'''

import argparse

from xsorb.cli.command import CLICommandBase


class CLICommand(CLICommandBase):
    """print how settings have been read from file,
    useful to check if parameters have been correctly set

Example:
 $ xsorb settings
    """

    @staticmethod
    def add_arguments(parser : argparse.ArgumentParser):
        parser.add_argument('-read_e_dft',
                            action='store_true',
                            help='if set, also try to read slab and molecule DFT energies')
        parser.add_argument('-read_e_ml',
                            action='store_true',
                            help='if set, also try to read slab and molecule ML energies')

    @staticmethod
    def run(args : argparse.Namespace):
        from xsorb.settings.settings import Settings # pylint: disable=import-outside-toplevel
        settings = Settings(read_energies_dft=args.read_e_dft,
                            read_energies_ml=args.read_e_ml)
        settings.print_summary()


    @staticmethod
    def bind_function(parser: argparse.ArgumentParser):
        parser.set_defaults(func=CLICommand.run)
