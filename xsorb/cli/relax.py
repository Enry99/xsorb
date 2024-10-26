'''
CLI parser for command: relax
'''

import argparse

from xsorb.cli.command import CLICommandBase, nonnegative_int, nonnegative_float


class CLICommand(CLICommandBase):
    """Select the configurations for the final relaxation, and launch the calculations

Examples:
 $ xsorb relax -n 10
 $ xsorb relax -t 0.1 --exclude 1 3 5
 $ xsorb relax -i 2 7 12 22
 $ xsorb relax -n 1 --by-site
 $ xsorb relax -n 1 --by-site --by-mol-idx
    """

    @staticmethod
    def add_arguments(parser : argparse.ArgumentParser):
        from xsorb.calculations.launchers import N_RELAX_DEFAULT
        relax_modes = parser.add_mutually_exclusive_group(required=True)
        relax_modes.add_argument('-n',
                                nargs='?',
                                type=nonnegative_int,
                                const=N_RELAX_DEFAULT,
                                default=None,
                                help='''Select the configurations for the final relaxation
                                by choosing the n lowest energy ones.
                                If --by-site is used, n configs are selected for each site''')
        relax_modes.add_argument('-t',
                                type=nonnegative_float,
                                help='''Select the configurations for the final relaxation
                                that have energy < min(energies) + t  (in eV)''')
        relax_modes.add_argument('-i',
                                nargs='+',
                                type=nonnegative_int,
                                help='''Manually select the configurations for the final relaxation
                                by their IDs''')
        parser.add_argument('--exclude',
                        nargs='+',
                        type=nonnegative_int,
                        help='''Exclude the configurations with the given IDs
                        from the final relaxation''')
        parser.add_argument('--by-site',
                        action='store_true',
                        default=False,
                        help='''Select the configurations for the final relaxation
                        on a per-site basis''')
        parser.add_argument('--by-mol-idx',
                        action='store_true',
                        default=False,
                        help='''Select the configurations for the final relaxation
                        on a per-molecule referenc atom basis''')
        parser.add_argument('--chem-phys',
                        action='store_true',
                        default=False,
                        help='''Select the configurations for the final relaxation
                        separating physisorbed and chemisorbed configurations.
                        The distinction is made based on a check based on covalent radii.
                        The multiplicative factor can be modified via the radius_scale_factor
                        key in the settings file. The default value is 1.1x''')
        parser.add_argument('--from',
                        type=str,
                        choices=['screening', 'ml_opt'],
                        default='screening',
                        dest='take_from',
                        help='''Choose the results on which perform the selection, whose
                        optimized structures will be used for the final relaxations, unless
                        --regenerate is specified. In that case, the structures will be re-generated
                        from scratch, and the -from argument will be used just to choose
                        which IDs should be included the''')
        parser.add_argument('--regenerate',
                        action='store_true',
                        help='''Re-generate the structures for the final relaxations
                        instead of reading the positions from the screening or ml_opt''')

    @staticmethod
    def run(args : argparse.Namespace):
        #validate args:
        if args.by_site or args.by_mol_idx or args.chem_phys and args.i:
            print('Error: cannot use --by-site, --by-mol-idx or --chem-phys with --i')
            import sys
            sys.exit(1)

        from xsorb.calculations.launchers import launch_final_relax
        launch_final_relax(n_configs=args.n,
                           threshold=args.t,
                           calc_ids=args.i,
                           excluded_calc_ids=args.exclude,
                           take_from=args.take_from,
                           relax_from_initial=args.regenerate,
                           by_site=args.by_site,
                           by_mol_idx=args.by_mol_idx,
                           separate_chem_phys=args.chem_phys,
                           )


    @staticmethod
    def bind_function(parser: argparse.ArgumentParser):
        parser.set_defaults(func=CLICommand.run)
