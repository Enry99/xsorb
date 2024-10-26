'''
Argument parser for xsorb command line interface.
'''

import argparse
from importlib import import_module

import xsorb
from xsorb.cli.command import CLICommandBase, CustomFormatter


HELP_MESSAGE="""

--------------------------------- Xsorb ---------------------------------

The program can automatically generate initial adsorption configurations combining
molecular rotations and adsorption sites, and optimize them by ab initio DFT calculations
(with VASP or Quantum Espresso) or by a machine learning model.

It employs a computationally efficient two(three)-step approach, with a first *screening* step
in which all the initial adsorption configurations are optimized with a larger convergence threshold,
and a second *relax* step, where only a selected subset of configurations is fully optimized.

The machine learning model can be employed:
i) as a pre-optmization tool, using the final positions as an initial guess for the screening procedue
ii) to substitute the ab initio screening
iii) as a single-run tool to fastly explore the configuration space.

--------------------------------- Paper ---------------------------------

Please cite the following paper if you use Xsorb in your research:
E. Pedretti, P. Restuccia, M.C. Righi, Comput Phys Commun, 291 (2023), Article 108827
https://doi.org/10.1016/j.cpc.2023.108827

------------------------------ Useful links -----------------------------

Official repository:    https://gitlab.com/triboteam/xsorbed
Latest updates:         https://github.com/Enry99/xsorb
Documentation:          https://gitlab.com/triboteam/xsorbed/-/wikis/home

"""


def build_xsorb_parser():
    '''
    Build the parser for the xsorb command line interface.
    '''

    # main parser
    parser = argparse.ArgumentParser(
        prog='xsorb',
        epilog=HELP_MESSAGE,
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
