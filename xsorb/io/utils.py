'''
General utility functions for the I/O module
'''

from ase.io import read, write

from xsorb import ase_custom
from xsorb.dft_codes.definitions import OUT_FILE_PATHS, LOG_FILE_PATHS, \
    SCF_NONCONVERGED_STRINGS, SCF_CONVERGED_STRINGS, OPTIMIZATION_COMPLETED_STRINGS


def ase_custom_read(filename, **kwargs):
    """
    Modified version of ase.io.read that returns an ase_custom.Atoms_custom object
    """
    atoms_or_atoms_list = read(filename=filename, **kwargs)
    if isinstance(atoms_or_atoms_list, list):
        return [ase_custom.Atoms_custom(at) for at in atoms_or_atoms_list]
    else:
        return ase_custom.Atoms_custom(atoms_or_atoms_list)


def overwrite_question(file_path : str) -> str:
    '''
    Asks the user if they want to overwrite the file at the specified path.

    Args:
    - file_path: str, path to the file

    Returns:
    - str, 'y' if the user wants to overwrite, 'n' otherwise,
        'yall' if the user wants to overwrite all files, 'nall' otherwise
    '''

    while True:
        answer = input(f'{file_path} already exists. Overwrite? '\
                       '("y" = yes to this one, "yall" = yes to all,'\
                         ' "n" = no to this one, "nall" = no to all): ')
        if answer.lower() in ['y', 'n', 'yall', 'nall']:
            return answer
        else:
            print('Value not recognized. Try again.')


def check_optimization_completed(program : str, calc_type : str, calc_id : int):
    '''
    Check if the given calculation is completed, reading the output file

    Args:
    - program: DFT program. Possible values: 'ESPRESSO' or 'VASP'
    - calc_type: 'SCREENING' or 'RELAX'
    - calc_id: numeric index of the calculation

    Returns:
    True or False
    '''

    filename = LOG_FILE_PATHS[calc_type][program].format(calc_id)
    searchfor = OPTIMIZATION_COMPLETED_STRINGS[program]

    with open(filename, 'r') as f:
        file_content = f.readlines()

    completed = False
    for line in file_content:
        if searchfor in line:
            completed = True
            break

    return completed


def scf_not_converged(program : str, calc_type : str, calc_id : int):
    '''
    Check if the given calculation has not reached SCF convergence, reading the output file

    Args:
    - program: DFT program. Possible values: 'ESPRESSO' or 'VASP'
    - calc_type: 'SCREENING' or 'RELAX'
    - calc_id: numeric index of the calculation

    Returns:
    True or False
    '''

    filename = OUT_FILE_PATHS[calc_type][program].format(calc_id)
    searchfor = SCF_NONCONVERGED_STRINGS[program]
    convergence_string = SCF_CONVERGED_STRINGS[program]

    with open(filename, 'r') as f:
        file_content = f.readlines()

    # we might encounter the situation where a first loop is not converged,
    # but the last one is, so we need to check all the lines:
    # the last one (conv or not conv) determines the status
    nonconv = False
    for line in file_content:
        if searchfor in line:
            nonconv = True
        elif convergence_string in line:
            nonconv = False

    return nonconv
