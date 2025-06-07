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
                        choices=['screening', 'relax', 'mlopt', 'all'],
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
        from xsorb.io.database import manual_update_calculations
        manual_update_calculations(args.calc_type,
                                   args.refresh,
                                   txt=args.txt)


    @staticmethod
    def bind_function(parser: argparse.ArgumentParser):
        parser.set_defaults(func=CLICommand.run)
