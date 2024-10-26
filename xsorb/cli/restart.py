'''
CLI parser for command: restart
'''

import argparse

from xsorb.cli.command import CLICommandBase


class CLICommand(CLICommandBase):
    """restart all (unfinished) screening or relax calculations

Example:
 $ xsorb restart screening
    """

    @staticmethod
    def add_arguments(parser : argparse.ArgumentParser):
        parser.add_argument('calc_type',
                        choices=['screening', 'relax'],
                        type=str,
                        help='''which type of calculation to restart''')

    @staticmethod
    def run(args : argparse.Namespace):
        from xsorb.io.jobs import restart_jobs
        restart_jobs(args.calc_type)


    @staticmethod
    def bind_function(parser: argparse.ArgumentParser):
        parser.set_defaults(func=CLICommand.run)
