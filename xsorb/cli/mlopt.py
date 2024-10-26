'''
CLI parser for command: mlopt
'''

import argparse

from xsorb.cli.command import CLICommandBase


class CLICommand(CLICommandBase):
    """generate the initial structures and launch the optimizations with a machine learning force field

    This can be used by itself, or as a pre-optimization step before the screening
    to have a better starting point for the DFT calculations, or even
    to substitute the screening step.

    """

    @staticmethod
    def add_arguments(parser : argparse.ArgumentParser):
        parser.add_argument('-save-images',
                        action='store_true',
                        help='''save an image with the adsorption sites,
                            and an image with the rotations of the molecule''')


    @staticmethod
    def run(args : argparse.Namespace):
        from xsorb.calculations.launchers import launch_ml_opt
        launch_ml_opt(save_images=args.save_images)


    @staticmethod
    def bind_function(parser: argparse.ArgumentParser):
        parser.set_defaults(func=CLICommand.run)
