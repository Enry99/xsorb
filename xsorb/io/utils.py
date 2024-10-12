from ase.io import read

from xsorb import ase_custom

def ase_custom_read(filename, **kwargs):
    return ase_custom.Atoms_custom(read(filename=filename, **kwargs))


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
