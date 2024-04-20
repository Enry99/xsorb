import argparse
from xsorbed.common_definitions import VERSION

def build_xsorb_parser():
    '''
    Build the argument list for the command-line parser

    Returns a ArgumentParser object
    '''
    
    parser = argparse.ArgumentParser(
        description='Select one of the commands below. The program will read further input informations from the file "settings.in"'
        )

    parser.add_argument('--v', '--version', action='version', version='%(prog)s '+VERSION)


    calc_group = parser.add_argument_group(title='Calculation options')
    #main flags
    main_calc = calc_group.add_mutually_exclusive_group()
    main_calc.add_argument('-g', action='store_true', help='generate all pwi(s) and a csv file with config labels, without submitting the jobs')
    main_calc.add_argument('-s', action='store_true', help='launch screening.')
    main_calc.add_argument('-r', action='store_true', help='launch final relaxations')
    main_calc.add_argument('-restart', choices=['s', 'r'], help='restart all (unfinished) screening or relax calculations')
    main_calc.add_argument('-regenerate_labels', action='store_true', help='Regenerate site_labels.csv if accidentally deleted.')
    #optional flags
    relax_options = calc_group.add_mutually_exclusive_group()
    relax_options.add_argument('--n', type=int, help='number of configurations (starting from the screening energy minimum) for the final relaxation. Default 5 (or 1 with --by-site)')
    relax_options.add_argument('--t', type=float, help='energy threshold above minimum (in eV) for choosing the configurations for the final relaxation')
    relax_options.add_argument('--i', nargs='+', type=int, help='select specific configurations (by label) for relaxation')
    calc_group.add_argument('--by-site', action='store_true', default=False, help='Select the most favorable configurations for full relax site-by-site site instead of overall')
    calc_group.add_argument('--exclude', nargs='+', type=int, help='exclude selected configurations from the final relaxation')
    calc_group.add_argument('--regenerate', action='store_true', default=False, help='re-generate structures for final relaxations instead of reading the positions from screening')
    calc_group.add_argument('--save-figs', action='store_true', help='plot also the figures (sites + molecule orientations) when using -g or -r')



    energy_group = parser.add_argument_group(title='Energies collection options')  
    #main flags
    main_energy = energy_group.add_mutually_exclusive_group()
    main_energy.add_argument('-e', action='store_true', help='save energies to file')
    energy_group.add_argument('--txt', action='store_true', help='write well formatted txt file instead of csv for -es or -er  (sorted by screen. en.)')
    main_energy.add_argument('-plot-energies-scr', action='store_true', help='save image with energies evolution during screening')
    main_energy.add_argument('-plot-energies', action='store_true', help='save image with energies evolution during final relax')


    visualization_group = parser.add_argument_group(title='Visualization options')
    #main flags
    main_visual = visualization_group.add_mutually_exclusive_group()
    main_visual.add_argument('-sites', action='store_true',  help='plot the identified high-symmetry sites sites with labels')
    main_visual.add_argument('-sites-all', action='store_true',  help='plot ALL the high-symmetry sites sites with labels')
    main_visual.add_argument('-screening-images', nargs='?', choices=['i', 'f'], const='f', help="save images of the screening configurations. Option: 'i'/'f' for initial or final positions. Default: f")
    main_visual.add_argument('-relax-images', nargs='?', choices=['i', 'f'], const='f', help="save images of the relaxations. Option: 'i'/'f' for initial or final positions. Default: f")
    #main_visual.add_argument('-render-image', nargs=2, help='select which image to render and the rotations list for the camera (syntax: [s/r][index] [rotations], e.g. r12 -10z,-80x)')
    main_visual.add_argument('-screening-animations', action='store_true', help='save animations of the screening in gif (default) or mp4 format (for povray)')
    main_visual.add_argument('-relax-animations', action='store_true', help='save animations of the full relaxations (screen.+final) in gif (default) or mp4 format (for povray)')
    main_visual.add_argument('-view', nargs=3, type=str, help='select which config to open with ase gui (syntax: [s/r] [index] [in/out], e.g. s 3 out).')
    main_visual.add_argument('-savefiles', nargs=3, help='save files in specific format (syntax: [format] [calc_type] [i/f], e.g. cif screening i or xyz relax f)')
    #optional flags
    visualization_group.add_argument('--povray', action='store_true', help='use povray to render images (if installed)')
    visualization_group.add_argument('--width-res', type=int, help='resolution width (in pixel) for povray')
    visualization_group.add_argument('--depth-cueing', nargs='?', type=float, const=1, help='Enable depth cueing. Optional parameter: intensity (>=0, default=1).')
    visualization_group.add_argument('--center-mol', action='store_true', help='translate the structure so that the molecule is centered in the image')
    visualization_group.add_argument('--cut-vacuum', action='store_true', help='cut most of the vacuum above molecule in the images')  
    visualization_group.add_argument('--rotation', type=str, help='Rotation for saving images, in ASE format, e.g. 10z,5x') 

    return parser


def validate_xsorb_args(args : argparse.Namespace, parser : argparse.ArgumentParser):
    '''
    Validate the arguments passed to the CLI, checking incompatibilities between the commands,
    and wrong format of the values passed to some options

    Args:
    - args: argparse.Namespace, returned by parse_args(), with the arguments read from command line
    - parser: argparse.ArgumentParser used to parse the arguments
    '''
    
    # Check presence of incompatible commands belonging to different groups
    args_dict = vars(args)
    calc_command = [ key for key in args_dict if args_dict[key] and (key=='s' or key=='g' or key=='r' or key=='restart') ]
    energy_command = [ key for key in args_dict if args_dict[key] and (key=='e' or key=='plot_energies_scr' or key=='plot_energies') ]
    visual_command = [ key for key in args_dict if args_dict[key] and (key=='sites' or key=='sites-all' or key=='screening-images'  
        or key=='relax-images' or key=='screening-animations' or key=='relax-animations' or key=='view' or key=='savefiles') ]

    if calc_command and energy_command:
        parser.error("argument -{1}: not allowed with argument -{0}".format(calc_command[0], energy_command[0]))
    elif calc_command and visual_command:
        parser.error("argument -{1}: not allowed with argument -{0}".format(calc_command[0], visual_command[0]))
    elif energy_command and visual_command:
        parser.error("argument -{1}: not allowed with argument -{0}".format(energy_command[0], visual_command[0]))
  

    # Check for all incompatible OPTIONAL arguments
    if args.save_figs and not args.g and not args.s:
        parser.error("--save-figs option can only be used with -g or -s")
    if args.n and not args.r:
        parser.error("--n option can only be used with -r")
    if args.t and not args.r:
        parser.error("--t option can only be used with -r")
    if args.i and not args.r:
        parser.error("--i option can only be used with -r")
    if args.by_site and not args.r:
        parser.error("--by-site option can only be used with -r")
    if args.by_site and args.i:
        parser.error("--by-site option cannot be used with --i")
    if args.exclude and not args.r:
        parser.error("--exclude option can only be used with -r")
    if args.regenerate and not args.r:
        parser.error("--regenerate option can only be used with -r")
    if args.povray and not args.screening_images and not args.relax_images and not args.screening_animations and not args.relax_animations:
        parser.error("--povray option can only be used with -screening-images or -relax-images or -screening-animations or -relax-animations")
    if args.width_res is not None and not args.povray:
        parser.error("--width-res can be specified only for povray rendering")
    if args.depth_cueing and not args.screening_images and not args.relax_images and not args.screening_animations and not args.relax_animations:
        parser.error("--depth-cueing option can only be used with -screening-images or -relax-images or -screening-animations or -relax-animations")
    if args.center_mol and not args.screening_images and not args.relax_images and not args.screening_animations and not args.relax_animations: # and not args.render_image:
        parser.error("--center-mol option can only be used with -screening-images or -relax-images or -screening-animations or -relax-animations")
    if args.rotation and not args.screening_images and not args.relax_images and not args.screening_animations and not args.relax_animations:
        parser.error("--rotation option can only be used with -screening-images or -relax-images or -screening-animations or -relax-animations")
    if args.cut_vacuum and not args.screening_images and not args.relax_images and not args.screening_animations and not args.relax_animations:
        parser.error("--cut-vacuum option can only be used with -screening-images or -relax-images or -screening-animations or -relax-animations")


    # Final check on values that are not already dealt with by argparse
    if(args.view):
        if args.view[0] != 's' and args.view[0] != 'r':
            parser.error('The first argument of -view must be either "s" or "r".')
        if(args.view[2] != 'in' and args.view[2] != 'out'):
            parser.error("The third argument of -view must be either 'in' or 'out'.")
        try:
            i = int(args.view[1:])
            if(i<0): raise ValueError("Index cannot be negative")
        except:
            parser.error('The second argument of -view must be a non-negative integrer. The syntax is: [s/r] [index] [in/out], e.g. s 3 out')


def cli_parse_xsorb():
    '''
    Parse the command-line arguments, and check that the inserted options are given in the correct way

    Returns the arguments in form of an argparse.Namespace
    '''

    parser = build_xsorb_parser()

    #parse the command line
    args = parser.parse_args()

    validate_xsorb_args(args=args, parser=parser)

    return args




def build_xfrag_parser():
    """
    Build the argument parser for the xfrag command-line interface.

    Returns:
        argparse.Namespace: The parsed command-line arguments.
    """

    parser = argparse.ArgumentParser(
        description='Select one of the commands below. The program will read further input informations from the files "fragments.in" and "settings.in"'
        )

    parser.add_argument('--v', '--version', action='version', version='%(prog)s '+VERSION)


    calc_group = parser.add_argument_group(title='Calculation options')
    #main flags
    main_calc = calc_group.add_mutually_exclusive_group()
    main_calc.add_argument('-generate-fragments', action='store_true', help='generate the input files of all framgents in fragments.json, wihtout launching calculations')
    main_calc.add_argument('-relax-fragments', action='store_true', help='generate the input files of all framgents listed in fragments.json and launch calculations')
    main_calc.add_argument('-g', action='store_true', help='generate all pwi(s) and a csv file with config labels, without submitting the jobs')
    main_calc.add_argument('-s', nargs='*', type=float, help='launch screening. Optional arguments: e_tot_conv_thr forc_conv_thr (default 5e-3 5e-2)')
    main_calc.add_argument('-r', action='store_true', help='launch final relaxations')
    main_calc.add_argument('-restart', choices=['s', 'r'], help='restart all (unfinished) screening or relax calculations')
    #optional flags
    relax_options = calc_group.add_mutually_exclusive_group()
    relax_options.add_argument('--n', type=int, help='number of configurations (starting from the screening energy minimum) for the final relaxation')
    relax_options.add_argument('--t', type=float, help='energy threshold above minimum (in eV) for choosing the configurations for the final relaxation')
    calc_group.add_argument('--by-site', action='store_true', default=False, help='Select the most favorable configurations for full relax site-by-site site instead of overall')


    energy_group = parser.add_argument_group(title='Energies collection options')  
    #main flags
    main_energy = energy_group.add_mutually_exclusive_group()
    main_energy.add_argument('-deltae', action='store_true', help='collect dissociation energies dE for all the (fully optimized) combinations (adsorbed molecule -> adsorbed fragments)')


    return parser


def validate_xfrag_args(args : argparse.Namespace, parser : argparse.ArgumentParser):
    '''
    Validate the command line arguments for the xfrag script.

    Args:
        args (argparse.Namespace): The parsed command line arguments.
        parser (argparse.ArgumentParser): The argument parser object.

    Returns:
        argparse.Namespace: The validated command line arguments.
    '''
    #Check presence of incompatible commands belonging to different groups
    args_dict = vars(args)
    calc_command = [ key for key in args_dict if args_dict[key] and (key == 'generate-fragments' or key == 'relax-fragments' or key=='s' or key=='g' or key=='r' or key=='restart') ]
    energy_command = [ key for key in args_dict if args_dict[key] and (key=='deltae') ]

    if calc_command and energy_command:
        parser.error(f"argument -{calc_command[0]}: not allowed with argument -{energy_command[0]}")  

    #Check for all incompatible OPTIONAL arguments
    if args.n and not args.r:
        parser.error("--n option can only be used with -r")
    if args.t and not args.r:
        parser.error("--t option can only be used with -r")
    if args.by_site and not args.r:
        parser.error("--by-site option can only be used with -r")

    #Final check on values that are not already dealt with by argparse
    if args.s is not None:
        if len(args.s) != 0 and len(args.s) != 2:
            parser.error('-s command requires either no argument (the default 5e-3 and 5e-2 will be used) or TWO arguments.')

    return args


def cli_parse_xfrag():
    '''
    Parse the command-line arguments, and check that the inserted options are given in the correct way

    Returns the arguments in form of an argparse.Namespace
    '''

    parser = build_xfrag_parser()

    #parse the command line
    args = parser.parse_args()

    validate_xfrag_args(args=args, parser=parser)

    return args