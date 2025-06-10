'''
CLI parser for command: gen
'''

import argparse

from xsorb.cli.command import CLICommandBase


class CLICommand(CLICommandBase):
    """Generate the initial structures, without launching any calculation.

The initial structures are generated based on the settings in the configuration file.
The structures are added to the structures.db database.

    """

    @staticmethod
    def add_arguments(parser : argparse.ArgumentParser):
        parser.add_argument('-save-images',
                        action='store_true',
                        help='''save an image with the adsorption sites,
                            and an image with the rotations of the molecule''')

    @staticmethod
    def run(args : argparse.Namespace):
        from xsorb.calculations.launchers import generate # pylint: disable=import-outside-toplevel
        generate(save_image=args.save_images)


    @staticmethod
    def bind_function(parser: argparse.ArgumentParser):
        parser.set_defaults(func=CLICommand.run)
