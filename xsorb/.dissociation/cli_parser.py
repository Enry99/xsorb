
# def build_xfrag_parser():
#     """
#     Build the argument parser for the xfrag command-line interface.

#     Returns:
#         argparse.Namespace: The parsed command-line arguments.
#     """

#     parser = argparse.ArgumentParser(
#         description='Select one of the commands below. The program will read further input informations from the files "fragments.in" and "settings.in"'
#         )

#     parser.add_argument('--v', '--version', action='version', version='%(prog)s '+VERSION)


#     calc_group = parser.add_argument_group(title='Calculation options')
#     #main flags
#     main_calc = calc_group.add_mutually_exclusive_group()
#     main_calc.add_argument('-generate-fragments', action='store_true', help='generate the input files of all framgents in fragments.json, wihtout launching calculations')
#     main_calc.add_argument('-relax-fragments', action='store_true', help='generate the input files of all framgents listed in fragments.json and launch calculations')
#     main_calc.add_argument('-g', action='store_true', help='generate all pwi(s) and a csv file with config labels, without submitting the jobs')
#     main_calc.add_argument('-preopt', action='store_true', help='Pre-optimize the structures before the screening using a machine learning force field.')
#     main_calc.add_argument('-s', action='store_true', help='launch screening.')
#     main_calc.add_argument('-r', action='store_true', help='launch final relaxations')
#     main_calc.add_argument('-restart', choices=['s', 'r'], help='restart all (unfinished) screening or relax calculations')
#     #optional flags
#     relax_options = calc_group.add_mutually_exclusive_group()
#     relax_options.add_argument('--n', type=int, help='number of configurations (starting from the screening energy minimum) for the final relaxation')
#     relax_options.add_argument('--t', type=float, help='energy threshold above minimum (in eV) for choosing the configurations for the final relaxation')
#     calc_group.add_argument('--by-site', action='store_true', default=False, help='Select the most favorable configurations for full relax site-by-site site instead of overall')
#     calc_group.add_argument('--from-preopt', action='store_true', help='Use the pre-optimized structures as starting configurations instead of the original ones')


#     energy_group = parser.add_argument_group(title='Energies collection options')
#     #main flags
#     main_energy = energy_group.add_mutually_exclusive_group()
#     main_energy.add_argument('-deltae', action='store_true', help='collect dissociation energies dE for all the (fully optimized) combinations (adsorbed molecule -> adsorbed fragments)')
#     main_energy.add_argument('-deltaeml', action='store_true', help='collect dissociation energies dE from the machine learning pre-optimization.')


#     return parser


# def validate_xfrag_args(args : argparse.Namespace, parser : argparse.ArgumentParser):
#     '''
#     Validate the command line arguments for the xfrag script.

#     Args:
#         args (argparse.Namespace): The parsed command line arguments.
#         parser (argparse.ArgumentParser): The argument parser object.

#     Returns:
#         argparse.Namespace: The validated command line arguments.
#     '''
#     #Check presence of incompatible commands belonging to different groups
#     args_dict = vars(args)
#     calc_command = [ key for key in args_dict if args_dict[key] and (key == 'generate-fragments' or key == 'relax-fragments' or key=='s' or key=='g' or key=='r' or key=='restart' or key=='preopt') ]
#     energy_command = [ key for key in args_dict if args_dict[key] and (key=='deltae') ]

#     if calc_command and energy_command:
#         parser.error(f"argument -{calc_command[0]}: not allowed with argument -{energy_command[0]}")

#     #Check for all incompatible OPTIONAL arguments
#     if args.n and not args.r:
#         parser.error("--n option can only be used with -r")
#     if args.t and not args.r:
#         parser.error("--t option can only be used with -r")
#     if args.by_site and not args.r:
#         parser.error("--by-site option can only be used with -r")

#     return args


# def cli_parse_xfrag():
#     '''
#     Parse the command-line arguments, and check that the inserted options are given in the correct way

#     Returns the arguments in form of an argparse.Namespace
#     '''

#     parser = build_xfrag_parser()

#     #parse the command line
#     args = parser.parse_args()

#     validate_xfrag_args(args=args, parser=parser)

#     return args