'''
Module for writing input files for the calculations
'''

from pathlib import Path

from xsorb.settings import Settings
from xsorb.structures import AdsorptionStructure
from xsorb.io.database import Database
from xsorb.io.utils import overwrite_question
from xsorb.dft_codes.definitions import IN_FILE_PATHS, OUT_FILE_PATHS
from xsorb.dft_codes.calculator import write_file_with_Calculator
from xsorb.dft_codes.override import override_dft_settings


def write_inputs(adsorption_structures : list[AdsorptionStructure],
                 settings : Settings,
                 calc_type : str | None = None,
                 override_settings : bool = True,
                 ask_before_overwrite : bool = True,
                 verbose : bool = True) -> list[int]:
    '''
    Writes the input files for all the adsorption configurations,
    updating the corresponding database(s).


    Args:
    - settings: Settings object, containing all the parameters
    - adsorption_structures: list AdsorptionStructure objects
    - calc_type: 'screening', 'relax' or 'ml_opt', or None (only generate input files)
    - override_settings: override some specifc settings (e.g. conv tresholds)
    - interactive: interactive mode: ask before overwriting files that are already present

    Returns:
    - written_indices: list of the indices of the configurations
        for which the input files have been written
    '''

    if verbose: print('Writing input files...') #pylint: disable=multiple-statements


    #DESIGN CHOICE: the database has ABSOLUTE PRIORITY over files.
    #When a new structure is created, it is ALWAYS written to the structures.db
    #If the user wants to play around a bit at the beginning to find the optimal parameters, they
    #can remove the structures.db file and start from scratch.
    #The other databases are updated only when not in generation mode.
    #The input files will be overwritten unless they are in their respective database.
    #In that case, the user will be asked if they want to overwrite the file or not.
    #If entry is in database and user wants to overwrite, check if the folder is present
    #and remove it before writing the new files.



    #when launching ml or screening from scratch, write to db_initial.
    #when launching screenining from preopt, read from db_ml, do not write to db_initial
    #when launching relax, read either from db_ml or db_screening.


    #Add structures to database, and obtain the calc_ids
    calc_ids = Database.add_structures(adsorption_structures)


    #Write the input files
    calc_type_for_writing = calc_type if calc_type is not None else 'screening'
    if override_settings:
        dftsettings = override_dft_settings(settings, program=settings.program, calc_type=calc_type)
    else:
        dftsettings = settings.dftprogram_settings_dict

    written_systems = []
    ANSWER_ALL = False
    for i, atoms in zip(calc_ids, adsorption_structures):

        file_label = f'{calc_type_for_writing.lower()}_{i}'  #e.g. screening_i or relax_i
        in_file_path = IN_FILE_PATHS[calc_type_for_writing][settings.program].format(i)
        out_file_path = OUT_FILE_PATHS[calc_type_for_writing][settings.program].format(i)
        file_dir = Path(in_file_path).parent

        if ask_before_overwrite and Path.exists(in_file_path) or Path.exists(out_file_path) \
            and not ANSWER_ALL:
            answer = overwrite_question(f'{in_file_path} or {out_file_path}')
            if answer in ('yall', 'nall'): ANSWER_ALL = True
            if answer in ('n', 'nall'): continue #skip if user does not want to overwrite

            #remove the directory
            Path.rm

        #initialize the Calculator and write input files
        write_file_with_Calculator(atoms=atoms,
                                   program=settings.program,
                                   dftsettings=dftsettings,
                                   label=file_label,
                                   directory=file_dir)

        written_systems.append({'calc_id': i,
                                'atoms': atoms,
                                'in_file_path': in_file_path,
                                'out_file_path': out_file_path})


    #if we are not in generation mode, update the databases.
    # if calc_id already present in db, but the file was rewritten,
    # rewrite also the row in the database.
    if calc_type is not None:







    if interactive: print('All input files written.')

    return written_indices


def write_inputs_ml(settings : Settings,
                 all_mol_on_slab_configs_ase : list,
                 calc_indices : list,
                 VERBOSE : bool = True):
    '''
    Writes the input files for all the adsorption configurations.


    Args:
    - settings: Settings object, containing the dft code parameters
    - all_mol_on_slab_configs_ase: list of ASE atoms for all the adsorption configurations

    Returns:
    - written_indices: list of the indices of the configurations for which the input files have been written
    '''

    if VERBOSE: print('Writing pre-optimization input files...')

    written_indices = []

    os.makedirs(preopt_outdir, exist_ok=True)

    for i, atoms in zip(calc_indices, all_mol_on_slab_configs_ase):

        os.makedirs(f'{preopt_outdir}/{i}', exist_ok=True)
        atoms.pbc = True #ensure that the periodic boundary conditions are set
        write(IN_FILE_PATHS['PREOPT']['ML'].format(i), atoms)
        written_indices.append(i)

    if VERBOSE: print('All input files written.')

    return written_indices