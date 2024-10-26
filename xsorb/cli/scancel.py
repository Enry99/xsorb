'''
CLI parser for command: scancel
'''

import argparse

from xsorb.cli.command import CLICommandBase


class CLICommand(CLICommandBase):
    """cancel all running jobs for this xsorb run (works only for Slurm scheduler)

    """

    @staticmethod
    def add_arguments(parser : argparse.ArgumentParser):
        return

    @staticmethod
    def run(args : argparse.Namespace):
        from xsorb.io.jobs import scancel
        scancel()


    @staticmethod
    def bind_function(parser: argparse.ArgumentParser):
        parser.set_defaults(func=CLICommand.run)
