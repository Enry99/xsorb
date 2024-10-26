'''
CLI parser for command: energyplot
'''

import argparse

from xsorb.cli.command import CLICommandBase


class CLICommand(CLICommandBase):
    """save image with energies evolution for a given calculation type

Example:
 $ xsorb energyplot screening
    """

    @staticmethod
    def add_arguments(parser : argparse.ArgumentParser):
        parser.add_argument('calc_type',
                        type=str,
                        choices=['screening', 'relax', 'ml_opt', 'initial'],
                        help='''type of calculation to render''')


    @staticmethod
    def run(args : argparse.Namespace):
        from xsorb.visualize.images import plot_energy_evolution
        plot_energy_evolution(args.calc_type)


    @staticmethod
    def bind_function(parser: argparse.ArgumentParser):
        parser.set_defaults(func=CLICommand.run)
