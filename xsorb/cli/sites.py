'''
CLI parser for command: sites
'''

import argparse

from xsorb.cli.command import CLICommandBase


class CLICommand(CLICommandBase):
    """plot the adsorption sites on the surface based on the current settings parameters.

Example:
 $ xsorb sites -all
    """

    @staticmethod
    def add_arguments(parser : argparse.ArgumentParser):
        parser.add_argument('-all',
                        action='store_true',
                        help='''plot all sites, regardless of the symmetry reduction
                        and the manual selection''')

    @staticmethod
    def run(args : argparse.Namespace):
        from xsorb.visualize.images import plot_adsorption_sites
        plot_adsorption_sites(all_sites=args.all)


    @staticmethod
    def bind_function(parser: argparse.ArgumentParser):
        parser.set_defaults(func=CLICommand.run)
