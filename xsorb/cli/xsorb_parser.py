'''
Argument parser for xsorb command line interface.
'''

import argparse
from importlib import import_module

import xsorb
from xsorb.cli.command import CLICommandBase, CustomFormatter


def build_xsorb_parser():
    '''
    Build the parser for the xsorb command line interface.
    '''

    # main parser
    parser = argparse.ArgumentParser(
        prog='xsorb',
        epilog=xsorb.__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        allow_abbrev=False)
    parser.add_argument('-v', '--version',action='version',version=f'%(prog)s-{xsorb.__version__}')
    parser.add_argument('-T', '--traceback',action='store_true',help='Print traceback on error')


    # subparsers
    subparsers = parser.add_subparsers(title='commands',dest='command')

    commands = [
        ('gen', 'xsorb.cli.gen'),
        ('screen', 'xsorb.cli.screen'),
        ('relax', 'xsorb.cli.relax'),
        ('mlopt', 'xsorb.cli.mlopt'),
        ('isolated', 'xsorb.cli.isolated'),
        ('restart', 'xsorb.cli.restart'),
        ('scancel', 'xsorb.cli.scancel'),
        ('sites', 'xsorb.cli.sites'),
        ('render', 'xsorb.cli.render'),
        ('view', 'xsorb.cli.view'),
        ('savefiles', 'xsorb.cli.savefiles'),
        ('update', 'xsorb.cli.update'),
        ('energyplot', 'xsorb.cli.energyplot')
    ]

    for command, module_name in commands:
        cmd : CLICommandBase = import_module(module_name).CLICommand

        subparser = subparsers.add_parser(
                    command,
                    formatter_class=CustomFormatter,
                    help=cmd.__doc__.split('\n')[0],
                    description=cmd.__doc__)
        cmd.add_arguments(subparser)
        cmd.bind_function(subparser)

    return parser
