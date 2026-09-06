'''
Module for setting up interactive calculators for different DFT codes based on user input.
'''

import os, sys

from ase import Atoms

def setup_calculator(code: str,
                     input_filename: str,
                     command: str|None,
                     pseudo_dir: str|None,
                     main_dir: str,
                     calc_id: int,
                     atoms : Atoms|None = None):
    '''
    Set up the calculator based on the specified code.

    Parameters
    ----------
    code : str
        The code to use for the calculation (e.g., 'espresso_unified', 'vasp_unified', 'ml').
    input_filename : str
        The filename of the input file containing the calculation settings.
    command : str
        The command to execute the calculation (e.g., 'mpirun -np 1 pw.x').
    pseudo_dir : str
        The directory containing the pseudopotentials for the calculation.
    main_dir : str
        The main directory containing the ml_calculator_loader.py (only needed for ML calculator).
    calc_id : int
        An identifier for the calculation, used for naming the socket.

    Returns
    -------
    calc : Calculator
        The initialized calculator object ready for use in the optimization.
    '''

    if code == 'espresso':

        # set up QE calculator
        from ase.calculators.espresso import Espresso, EspressoProfile
        from xsorb.ase_custom import espresso # patched version of espresso calculator
        from xsorb.dft_codes.input_settings import build_espresso_settings_dict

        dftsettings = build_espresso_settings_dict(input_filename)
        pseudopotentials = dftsettings.pop('pseudopotentials')

        profile = EspressoProfile(
            command=command,
            pseudo_dir=pseudo_dir
        )

        calc = Espresso(
            profile=profile,
            directory='calc',
            pseudopotentials=pseudopotentials,
            kpts=dftsettings["kpts"],
            koffset=dftsettings["koffset"],
            input_data=dftsettings,
        ).socketio(unixsocket=f'espresso-xsorb-{calc_id}')

    ###########################################################################
    elif code == "vasp":

        # set up VASP calculator
        from vasp_interactive import VaspInteractive

        os.environ["VASP_PP_PATH"] = pseudo_dir

        calc = VaspInteractive(command=command,
                            atoms=atoms,
                            directory='.')
        calc.read_incar('INCAR')
        calc.read_kpoints('KPOINTS')
        calc.read_potcar('POTCAR')

    ###########################################################################
    elif code == 'ml':
        try:
            sys.path.append(os.path.abspath(main_dir))
            from ml_calculator_loader import NNloader #runtime import
        except ImportError as exc:
            raise ImportError("ml_calculator_loader.py not found.") from exc

        try:
            calc = NNloader()
        except Exception as exc:
            raise RuntimeError("Error loading the ML calculator.") from exc

    else:
        raise ValueError(f"Unsupported code: {code}")

    return calc
