'''
General utility functions for the I/O module
'''

from ase.io import read, write

from xsorb import ase_custom


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


#TODO: implement this function
# def saveas(calc_type : str, i_or_f : str, saveas_format : str):
#     '''
#     Save all the configurations in a different format, e.g. xyz or cif.

#     Args:
#     - calc_type: 'screening' or 'relax'
#     - i_or_f: initial or final coordinates of the relaxation (both for screening and full relax)
#     - saveas_format: file format, e.g. xyz
#     '''

#     settings = Settings(read_energies=False)

#     if i_or_f == 'i':
#         FILE_PATHS = IN_FILE_PATHS
#     elif i_or_f == 'f':
#         FILE_PATHS = OUT_FILE_PATHS
#     else:
#         raise RuntimeError(f"Wrong arguments: passed '{calc_type} {i_or_f}', expected 'screening i/f' or 'relax i/f'")
#     if calc_type != 'screening' and calc_type != 'relax':
#         raise RuntimeError(f"Wrong argument: passed '{calc_type}', expected 'screening' or 'relax'")

#     folder = f"{saveas_format}/{calc_type}"

#     print(f"Saving files to {folder}...")
#     os.makedirs(folder, exist_ok=True)

#     indices = _get_configurations_numbers()
#     for i in indices:
#         if os.path.isfile(FILE_PATHS[calc_type.upper()][settings.program].format(i)):
#             atoms = read(FILE_PATHS[calc_type.upper()][settings.program].format(i))
#             if(saveas_format == 'xyz'):
#                 ase_custom.write_xyz_custom(f'{folder}/{calc_type}_{i}.{saveas_format}', atoms)
#             else:
#                 write(f'{folder}/{calc_type}_{i}.{saveas_format}', atoms)

#     print("All files saved.")