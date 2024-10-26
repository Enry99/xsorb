'''
CLI parser for command: screen
'''

import argparse

from xsorb.cli.command import CLICommandBase


class CLICommand(CLICommandBase):
    """generate the initial structures and launch the screening

    The structures are generated based on the input parameters and the adsorption sites.
    For the screening, higher convergence thresholds are used to speed up the calculations.
    The value(s) for the convergence thresholds can be set in the input file, otherwise
    the default values are used (see the documentation).
    When launching the screening, the screening.db database is created in the current directory.

    """

    @staticmethod
    def add_arguments(parser : argparse.ArgumentParser):
        parser.add_argument('-from-ml',
                        action='store_true',
                        help='''Use the final structures from machine learning optimization
                        as starting points for the screening''')
        parser.add_argument('-save-images',
                        action='store_true',
                        help='''save an image with the adsorption sites,
                            and an image with the rotations of the molecule''')

    @staticmethod
    def run(args : argparse.Namespace):
        from xsorb.calculations.launchers import launch_screening
        launch_screening(from_ml_opt=args.from_ml, save_image=args.save_images)


    @staticmethod
    def bind_function(parser: argparse.ArgumentParser):
        parser.set_defaults(func=CLICommand.run)
