#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Enrico Pedretti

'''
Module for writing input files for the calculations
'''

from __future__ import annotations
import shutil
from pathlib import Path
from dataclasses import dataclass

from ase import Atoms

from xsorb.structures.generation import AdsorptionStructure
from xsorb.structures.utils import set_fixed_slab_constraints
from xsorb.io.settings import Settings
from xsorb.io.database import Database
from xsorb.io.utils import overwrite_question, write, write_xyz_custom
from xsorb.dft_codes.definitions import IN_FILE_PATHS, OUT_FILE_PATHS, LOG_FILE_PATHS
from xsorb.dft_codes.calculator import write_file_with_calculator
from xsorb.dft_codes.override import override_dft_settings


@dataclass
class WrittenSystem:
    '''
    Small dataclass to store info about the written systems
    '''
    calc_id: int | str #index or 'slab'/'mol'
    adsorption_structure: AdsorptionStructure | None
    in_file_path: str
    out_file_path: str
    log_file_path: str
    job_id: int | None = None
    adsite_z: float | None = None
    mol_ref_idx: int | None = None


def write_inputs(*,adsorption_structures : list[AdsorptionStructure],
                 settings : Settings,
                 calc_type : str | None = None,
                 calc_ids : list[int] | None = None,
                 override_settings : bool = True,
                 ask_before_overwrite : bool = True,
                 verbose : bool = True) -> list[WrittenSystem]:
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
    - written_systems: list of WrittenSystem objects, containing the calc_id,
        the AdsorptionStructure object, the path to the input, output and log files.
    '''

    program = settings.program if calc_type is not 'ml_opt' else 'ml'

    if verbose: print('Writing input files...') #pylint: disable=multiple-statements

    #Add structures to database, and obtain the calc_ids
    if calc_ids is None:
        calc_ids = Database.add_structures(adsorption_structures)


    #Write the input files
    calc_type_for_writing = calc_type if calc_type is not None else 'screening'
    if override_settings and program != 'ml' and calc_type is not None:
        dftsettings = override_dft_settings(settings.dftprogram_settings_dict,
                                            program=program,
                                            calc_type=calc_type)
    else:
        dftsettings = settings.dftprogram_settings_dict

    written_systems = []
    answer_all = False #pylint: disable=invalid-name
    for i, ads_structure in zip(calc_ids, adsorption_structures):

        #possibly apply constraints to slab in case of ml_opt
        if calc_type is 'ml_opt' and settings.structure.constraints.fix_slab_ml_opt:
            set_fixed_slab_constraints(ads_structure.atoms, ads_structure.slab_indices)

        file_label = f'{calc_type_for_writing.lower()}_{i}'  #e.g. screening_i or relax_i
        in_file_path = IN_FILE_PATHS[calc_type_for_writing][program].format(i)
        out_file_path = OUT_FILE_PATHS[calc_type_for_writing][program].format(i)
        log_file_path = LOG_FILE_PATHS[calc_type_for_writing][program].format(i)
        file_dir = Path(in_file_path).parent.as_posix()

        if ask_before_overwrite and Path(in_file_path).exists() or Path(out_file_path).exists() \
            and not answer_all:
            answer = overwrite_question(f'{in_file_path} or {out_file_path}')
            if answer in ('yall', 'nall'): answer_all = True #pylint: disable=multiple-statements,invalid-name
            if answer in ('n', 'nall'): continue  #pylint: disable=multiple-statements

            #remove the directory and all its content
            shutil.rmtree(file_dir)

        #initialize the Calculator and write input files
        write_file_with_calculator(atoms=ads_structure.atoms,
                                   program=program,
                                   dftsettings=dftsettings,
                                   label=file_label,
                                   directory=file_dir)

        written_systems.append(WrittenSystem(calc_id=i,
                                             adsorption_structure=ads_structure,
                                             in_file_path=in_file_path,
                                             out_file_path=out_file_path,
                                             log_file_path=log_file_path))


    #if we are not in generation mode, update the databases.
    if calc_type is not None:
        total_e_slab_mol=settings.input.total_e_slab_mol if program != 'ml' else None
        Database.add_calculations(systems=written_systems,
                                  program=program,
                                  mult=settings.structure.molecule.radius_scale_factor,
                                  total_e_slab_mol=total_e_slab_mol,
                                  calc_type=calc_type)

    if verbose: print('All input files written.') #pylint: disable=multiple-statements

    return written_systems


def write_slab_mol_inputs(*,slab : Atoms | None,
                          molecule : Atoms | None,
                          settings : Settings,
                          ml : bool,
                          ask_before_overwrite : bool = True,
                          verbose : bool = True) -> list[WrittenSystem]:
    '''
    Writes the input files for slab and/or molecule.

    Args:
    - slab: Atoms object for the slab
    - molecule: Atoms object for the molecule
    - settings: Settings object, containing all the parameters
    - ml: True if the input files are for the machine learning model
    - ask_before_overwrite: interactive mode: ask before overwriting files that are already present

    Returns:
    - written_systems: list of WrittenSystem objects, containing the calc_id,
        the Atoms object, the path to the input, output and log files.
    '''

    program = 'ml' if ml else settings.program
    dftsettings = settings.dftprogram_settings_dict

    structures, written_systems = [], []
    if slab is not None:
        structures.append(slab)
        written_systems.append(WrittenSystem(calc_id='slab',
                                             adsorption_structure=None,
                                             in_file_path=IN_FILE_PATHS['slab'][program],
                                             out_file_path=OUT_FILE_PATHS['slab'][program],
                                             log_file_path=LOG_FILE_PATHS['slab'][program]))

    if molecule is not None:
        structures.append(molecule)
        written_systems.append(WrittenSystem(calc_id='mol',
                                             adsorption_structure=None,
                                             in_file_path=IN_FILE_PATHS['mol'][program],
                                             out_file_path=OUT_FILE_PATHS['mol'][program],
                                             log_file_path=LOG_FILE_PATHS['mol'][program]))


    if verbose: print('Writing input files...') #pylint: disable=multiple-statements

    #Write the input files
    written_systems = []
    answer_all = False #pylint: disable=invalid-name
    for atoms, system in zip(structures, written_systems):

        #possibly apply constraints to slab in case of ml_opt
        if ml and system.calc_id is 'slab' and settings.structure.constraints.fix_slab_ml_opt:
            set_fixed_slab_constraints(atoms)

        file_label : str = system.calc_id
        in_file_path = system.in_file_path
        out_file_path = system.out_file_path
        file_dir = Path(in_file_path).parent.as_posix()

        if ask_before_overwrite and Path(in_file_path).exists() or Path(out_file_path).exists() \
            and not answer_all:
            answer = overwrite_question(f'{in_file_path} or {out_file_path}')
            if answer in ('yall', 'nall'): answer_all = True #pylint: disable=multiple-statements,invalid-name
            if answer in ('n', 'nall'): continue  #pylint: disable=multiple-statements

            #remove the directory and all its content
            shutil.rmtree(file_dir)

        #initialize the Calculator and write input files
        write_file_with_calculator(atoms=atoms,
                                   program=program,
                                   dftsettings=dftsettings,
                                   label=file_label,
                                   directory=file_dir)

    if verbose: print('All input files written.') #pylint: disable=multiple-statements

    return written_systems


def saveas(calc_type : str, saveas_format : str):
    '''
    Save all the configurations in a different format, e.g. xyz or cif.
    If initial, the configurations are read from structures.db, and correspond
    to the initial generated ones. Otherwise, the final configurations for
    the given mode are written.

    Args:
    - calc_type: 'initial','screening','relax','ml_opt'
    - saveas_format: file format, e.g. xyz
    '''

    if calc_type not in ('initial','screening', 'relax', 'ml_opt'):
        raise RuntimeError(f"Wrong '{calc_type}', expected 'screening', 'relax' or 'ml_opt'")

    folder = Path(f"{calc_type}/{saveas_format}")

    print(f"Saving files to {folder}...")
    folder.mkdir(exist_ok=True, parents=True)

    if calc_type == 'initial':
        rows = Database.get_structures()
    else:
        rows = Database.get_calculations(calc_type)

    for row in rows:
        if saveas_format == 'xyz':
            write_xyz_custom(folder / f'{calc_type}_{row.calc_id}.{saveas_format}', row.atoms)
        else:
            write(folder / f'{calc_type}_{row.calc_id}.{saveas_format}', row.atoms)

    print("All files saved.")
