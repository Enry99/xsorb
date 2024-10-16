'''
General utility functions for the I/O module
'''

from ase.io import read, write

from xsorb.ase_custom import write_xyz_custom, Atoms_custom


def ase_custom_read(filename, **kwargs):
    """
    Modified version of ase.io.read that returns an ase_custom.Atoms_custom object
    """
    atoms_or_atoms_list = read(filename=filename, **kwargs)
    if isinstance(atoms_or_atoms_list, list):
        return [Atoms_custom(at) for at in atoms_or_atoms_list]
    else:
        return Atoms_custom(atoms_or_atoms_list)


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


def continue_even_if_not_all_completed_question() -> bool:
    '''
    Asks the user if they want to continue even if not all calculations are completed.

    Returns:
    - bool, True if the user wants to continue, False otherwise
    '''

    while True:
        answer = input('Not all calculations are completed. '\
                       'Continue anyway with those that are present? (y/n): ')
        if answer.lower() in ['y', 'n']:
            return answer.lower() == 'y'
        else:
            print('Value not recognized. Try again.')
