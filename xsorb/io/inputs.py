'''
Module for writing input files for the calculations
'''

import shutil
from pathlib import Path

from xsorb.settings import Settings
from xsorb.structures import AdsorptionStructure, set_fixed_slab_constraints
from xsorb.io.database import Database
from xsorb.io.utils import overwrite_question
from xsorb.dft_codes.definitions import IN_FILE_PATHS, OUT_FILE_PATHS
from xsorb.dft_codes.calculator import write_file_with_Calculator
from xsorb.dft_codes.override import override_dft_settings


def write_inputs(adsorption_structures : list[AdsorptionStructure],
                 settings : Settings,
                 calc_type : str | None = None,
                 calc_ids : list[int] | None = None,
                 override_settings : bool = True,
                 ask_before_overwrite : bool = True,
                 verbose : bool = True) -> list[dict]:
    '''
    Writes the input files for all the adsorption configurations,
    updating the corresponding database(s).


    Args:
    - adsorption_structures: list of AdsorptionStructure objects
    - settings: Settings object, containing all the parameters
    - calc_type: 'screening', 'relax' or 'ml_opt', or None (only generate input files)
    - calc_ids: list of the calculation IDs. If None, the IDs are automatically assigned
        They will be None when called from the generation mode (from scratch)
    - override_settings: override some specifc settings (e.g. conv tresholds)
    - interactive: interactive mode: ask before overwriting files that are already present

    Returns:
    - written_systems: list of dictionaries, containing the calc_id,
        the AdsorptionStructure object, the path to the input file and the path to the output file
    '''

    if verbose: print('Writing input files...') #pylint: disable=multiple-statements

    #Add structures to database, and obtain the calc_ids
    if calc_ids is None:
        calc_ids = Database.add_structures(adsorption_structures)


    #Write the input files
    calc_type_for_writing = calc_type if calc_type is not None else 'screening'
    if override_settings:
        dftsettings = override_dft_settings(settings, program=settings.program, calc_type=calc_type)
    else:
        dftsettings = settings.dftprogram_settings_dict

    written_systems = []
    ANSWER_ALL = False #pylint: disable=invalid-name
    for i, ads_structure in zip(calc_ids, adsorption_structures):

        #possibly apply constraints to slab in case of ml_opt
        if calc_type is 'ml_opt' and settings.structure.constraints.fix_slab_preopt:
            set_fixed_slab_constraints(ads_structure.atoms, ads_structure.slab_indices)

        file_label = f'{calc_type_for_writing.lower()}_{i}'  #e.g. screening_i or relax_i
        in_file_path = IN_FILE_PATHS[calc_type_for_writing][settings.program].format(i)
        out_file_path = OUT_FILE_PATHS[calc_type_for_writing][settings.program].format(i)
        file_dir = Path(in_file_path).parent

        if ask_before_overwrite and Path.exists(in_file_path) or Path.exists(out_file_path) \
            and not ANSWER_ALL:
            answer = overwrite_question(f'{in_file_path} or {out_file_path}')
            if answer in ('yall', 'nall'): ANSWER_ALL = True #pylint: disable=multiple-statements,invalid-name
            if answer in ('n', 'nall'): continue  #pylint: disable=multiple-statements

            #remove the directory and all its content
            shutil.rmtree(file_dir)

        #initialize the Calculator and write input files
        write_file_with_Calculator(atoms=ads_structure.atoms,
                                   program=settings.program,
                                   dftsettings=dftsettings,
                                   label=file_label,
                                   directory=file_dir)

        written_systems.append({'calc_id': i,
                                'adsorption_structure': ads_structure,
                                'in_file_path': in_file_path,
                                'out_file_path': out_file_path})


    #if we are not in generation mode, update the databases.
    if calc_type is not None:
        Database.add_calculations(written_systems, calc_type)

    if verbose: print('All input files written.') #pylint: disable=multiple-statements

    return written_systems



def write_slab_mol_ml(slab, mol):
    os.makedirs(preopt_outdir, exist_ok=True)
    os.makedirs(preopt_outdir+'/slab', exist_ok=True)
    os.makedirs(preopt_outdir+'/mol', exist_ok=True)

    write(preopt_outdir+'/slab/slab.xyz', slab)
    write(preopt_outdir+'/mol/mol.xyz', mol)
