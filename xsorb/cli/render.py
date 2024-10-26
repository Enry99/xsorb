'''
CLI parser for command: render
'''

import argparse

from xsorb.cli.command import CLICommandBase, nonnegative_int, nonnegative_float


class CLICommand(CLICommandBase):
    """render images of all the configurations for a given calculation type

Example:
 $ xsorb render screening -pov
 $ xsorb render relax 12
    """

    @staticmethod
    def add_arguments(parser : argparse.ArgumentParser):
        parser.add_argument('calc_type',
                        type=str,
                        choices=['screening', 'relax', 'ml_opt', 'initial'],
                        help='''type of calculation to render''')
        parser.add_argument('calc_id',
                        nargs='?',
                        type=nonnegative_int,
                        help='''ID of the calculation to render''')
        parser.add_argument('-pov','--povray',
                        action='store_true',
                        default=False,
                        help='Use povray for rendering (much better quality).')
        parser.add_argument('-s','--supercell',
                        nargs = 3,
                        type=nonnegative_int,
                        metavar=('nx', 'ny', 'nz'),
                        help="Replicate the cell nx ny nz times along the three cell vectors.")
        parser.add_argument('-wr', '--wrap',
                        action='store_true',
                        default=False,
                        help='Wrap atoms according to pbc.')
        parser.add_argument('-cm', '--center-mol',
                        action='store_true',
                        help='translate the structure so that the molecule is centered in the image')
        parser.add_argument('-r', '--rotation',
                        type=str,
                        help='Rotation for saving images, in ASE format, e.g. 10z,5x')
        parser.add_argument('-dc','--depth-cueing',
                        nargs='?',
                        type=nonnegative_float,
                        const=1.0,
                        help='Enable depth cueing. Optional parameter: intensity (>0, default=1).')
        parser.add_argument('-rc','--range-cut',
                        nargs=2,
                        type=float,
                        metavar=('zmin', 'zmax'),
                        help='range to be displayed in the z direction.')
        parser.add_argument('-cv','--cut-vacuum',
                        action='store_true',
                        help='Cut vacuum above and below the slab (avoid white empty region).')
        parser.add_argument('-cc', '--colorcode',
                        type=str,
                        choices=['forces', 'magmoms', 'coordnum'],
                        help='''Color atoms according to a property
                        (e.g. forces, magnetic moments, coordination number).''')
        parser.add_argument('--ccrange',
                        nargs=2,
                        type=float,
                        metavar=('min', 'max'),
                        help='range for the colorcode.')
        parser.add_argument('-arr', '--arrows',
                        type=str,
                        choices=['forces', 'magmoms'],
                        help='''Draw arrows representing the vectors,
                        with lenghth proportional to the magnitude.''')
        parser.add_argument('-nobd','--nobonds',
                        action='store_true',
                        default=False, help='Do not show bonds' )
        parser.add_argument('-w', '--width-resolution',
                        type=nonnegative_int,
                        default=700,
                        help='Horizontal resolution in pixels.')
        parser.add_argument('-m','--movie',
                        action='store_true',
                        default=False,
                        help='Create movie from the frames.')
        parser.add_argument('-f', '--framerate',
                        type=nonnegative_int,
                        default=10,
                        help='Framerate of the movie (frames per second). Default = 10.')

    @staticmethod
    def run(args : argparse.Namespace):
        args_dict = vars(args)
        calc_type = args_dict.pop('calc_type')
        calc_id = args_dict.pop('calc_id')
        movie = args_dict.pop('movie')

        from xsorb.visualize.images import plot_images
        plot_images(calc_type=calc_type, calc_id=calc_id, movie=movie, **args_dict)


    @staticmethod
    def bind_function(parser: argparse.ArgumentParser):
        parser.set_defaults(func=CLICommand.run)
