'''
CLI parser for command: isolated
'''

import argparse

from xsorb.cli.command import CLICommandBase


class CLICommand(CLICommandBase):
    """launch the calculations for the isolated slab, molecule or both.

    Can be used for DFT or ML calculations.

    Example:
    xsorb isolated -slab -mol

    """

    @staticmethod
    def add_arguments(parser : argparse.ArgumentParser):
        parser.add_argument('-ml',
                        action='store_true',
                        help='''launch the calculation with machine learning potential''')
        parser.add_argument('-slab',
                        action='store_true',
                        help='''launch the calculation for the isolated slab''')
        parser.add_argument('-molecule',
                        action='store_true',
                        help='''launch the calculation for the isolated molecule''')
        parser.add_argument('-samecell',
                        action='store_true',
                        help='''use the same cell of the slab for the molecule.
                        If not set and the molecule has no cell, an orthorombic cell
                        with 10 Angstrom of vacuum in each direction will be used ''')


    @staticmethod
    def run(args : argparse.Namespace):
        from xsorb.calculations.launchers import launch_isolated_slab_and_molecule
        launch_isolated_slab_and_molecule(ml=args.ml,
                                          launch_slab=args.slab,
                                          launch_mol=args.molecule,
                                          samecell=args.samecell)


    @staticmethod
    def bind_function(parser: argparse.ArgumentParser):
        parser.set_defaults(func=CLICommand.run)
