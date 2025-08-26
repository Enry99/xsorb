'''
CLI parser for command: isolated
'''

import argparse

from xsorb.cli.command import CLICommandBase


class CLICommand(CLICommandBase):
    """launch the calculations for the isolated slab, molecule or both.

    Can be used for DFT (default) or ML (with option -ml) calculations.

    Note that as of now this is used just for obtaining slab/molecule reference energies,
    while the structures are always taken from the provided files.

    Example:
    $ xsorb isolated
    $ xsorb isolated -ml
    $ xsorb isolated molecule -samecell

    """

    @staticmethod
    def add_arguments(parser : argparse.ArgumentParser):
        parser.add_argument('which',
                        type=str,
                        nargs='?',
                        choices=['slab', 'molecule', 'both'],
                        default='both',
                        help='''launch the calculation for the isolated slab, molecule or both''')
        parser.add_argument('-ml',
                        action='store_true',
                        help='''use machine learning potential (default is DFT)''')
        parser.add_argument('-samecell',
                        action='store_true',
                        help='''use the same cell of the slab for the molecule.
                        If not set and the molecule has no cell, an orthorombic cell
                        with 10 Angstrom of vacuum in each direction will be used ''')
        parser.add_argument('-use-constraints',
                            action='store_true',
                            help='''include the constraints defined in settings file''')


    @staticmethod
    def run(args : argparse.Namespace):
        from xsorb.calculations.launchers import launch_isolated_slab_and_molecule # pylint: disable=import-outside-toplevel
        launch_isolated_slab_and_molecule(ml=args.ml,
                                          launch_slab=args.which in ['slab', 'both'],
                                          launch_mol=args.which in ['molecule', 'both'],
                                          samecell=args.samecell,
                                          use_constraints=args.use_constraints)


    @staticmethod
    def bind_function(parser: argparse.ArgumentParser):
        parser.set_defaults(func=CLICommand.run)
