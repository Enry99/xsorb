'''
CLI parser for command: update
'''

import argparse

from xsorb.cli.command import CLICommandBase


class CLICommand(CLICommandBase):
    """update the database with the latest results fom the calculations

    Refresh to force the database to update all the calculations,
    re-reading all the output files.

    If you changed the radius_scale_factor in the settings,
    you should refresh the database to update the bonding status.

    """

    @staticmethod
    def add_arguments(parser : argparse.ArgumentParser):
        parser.add_argument('calc_type',
                        type=str,
                        choices=['screening', 'relax', 'ml_opt', 'all'],
                        help='''type of calculation to plot''')
        parser.add_argument('-refresh',
                        action='store_true',
                        default=False,
                        help='''refresh the database''')
        parser.add_argument('--txt',
                        action='store_true',
                        default=False,
                        help='''output to txt instead of csv''')


    @staticmethod
    def run(args : argparse.Namespace):
        from xsorb.io.database import Database

        if args.refresh:
            from xsorb.io.settings import Settings
            settings = Settings(read_energies=True)
            total_e_slab_mol = settings.total_e_slab_mol
            mult=settings.structure.molecule.radius_scale_factor
        else:
            mult = None
            total_e_slab_mol = None
        Database.update_calculations(args.calc_type,
                                     args.refresh,
                                     total_e_slab_mol=total_e_slab_mol,
                                     mult=mult,
                                     txt=args.txt,
                                     verbose=True)


    @staticmethod
    def bind_function(parser: argparse.ArgumentParser):
        parser.set_defaults(func=CLICommand.run)
