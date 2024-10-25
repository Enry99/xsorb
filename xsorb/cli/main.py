'''
Module for parsing command line arguments.
'''

import argparse
from importlib import import_module

#import xsorb

#TODO: in the whole code, replace strings \ with multiline strings (docstrings)

parser = argparse.ArgumentParser(
    prog='xsorb',
    #description=xsorb.__doc__,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    allow_abbrev=False
)

#parser.add_argument('-v', '--version', action='version', version=f'%(prog)s-{xsorb.__version__}')

subparsers = parser.add_subparsers(title='Sub-commands',dest='command')

commands = [
    ('gen', 'gen'),
    # 'screen',
    # 'relax',
    # 'mlopt',
    # 'restart',
    # 'scancel',
    #'sites',
    #'regenerate_labels',
]

for command, module_name in commands:
    cmd = import_module(module_name).CLICommand
    subparser = subparsers.add_parser(
                command,
                formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                help=short,
                description=long)
    cmd.add_arguments(subparser)
    cmd.bind_function(subparser)




# def sites(args):
#     return
#     from xsorb.visualize import plot_adsorption_sites
#     plot_adsorption_sites(all_sites=args.all)

# sites_parser = subparsers.add_parser('sites',
#     help='plot the adsorption sites on the surface identified the current settings parameters')
# sites_parser.set_defaults(func=sites)
# sites_parser.add_argument('-a', '--all',
#                           action='store_true',
#                           help='show the plot')




screen_parser = subparsers.add_parser('screen', help='generate the initial structures and launch the screening')
relax_parser = subparsers.add_parser('relax', help='select the lowest energy configurations from screening or ml optimization and launch the final relaxations')
mlopt_parser = subparsers.add_parser('mlopt', help='pre-optimize the structures before the screening using a machine learning force field')
restart_parser = subparsers.add_parser('restart', help='restart all (unfinished) screening or relax calculations')
scancel_parser = subparsers.add_parser('scancel', help='cancel all running jobs for this xsorb run (works only for Slurm scheduler)')
#subparsers.add_parser('regenerate_labels')

#subparsers = parser.add_subparsers(title='visualization commands')





args = parser.parse_args("gen".split())
args.func(args)

#print(args)



# def hyphenated(string):
#     return '-'.join([word[:4] for word in string.casefold().split()])

# parser = argparse.ArgumentParser()
# _ = parser.add_argument('short_title', type=hyphenated)


# parser = argparse.ArgumentParser()
# parser.add_argument('--foo', metavar='YYY')
# parser.add_argument('bar', metavar='XXX')
# parser.parse_args('X --foo Y'.split())
# parser.print_help()
# #usage:  [-h] [--foo YYY] XXX