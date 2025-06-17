'''
CLI parser for command: cleanup
'''

import argparse

from xsorb.cli.command import CLICommandBase


class CLICommand(CLICommandBase):
    """Remove all files for an Xsorb run.
    If `calc_only` is set, only calculation files are removed,
    keeping the databases (and adsites file).
    """

    @staticmethod
    def add_arguments(parser : argparse.ArgumentParser):
        parser.add_argument(
            '--calc_only',
            action='store_true',
            help='Clean only calculation files, keeping the databases and adsites file.'
        )

    @staticmethod
    def run(args : argparse.Namespace):
        from xsorb.io.cleanup import cleanup_xsorb_run  # pylint: disable=import-outside-toplevel
        cleanup_xsorb_run(calc_only=args.calc_only)


    @staticmethod
    def bind_function(parser: argparse.ArgumentParser):
        parser.set_defaults(func=CLICommand.run)
