#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue 28 Feb 2023
@author: Enrico Pedretti

Function definitions to read from pwo and launch scripts

"""

import os, shutil
import numpy as np
import pandas as pd
from ase.io import read, write
from dftcode_specific import edit_files_for_restart, UNITS_TO_EV_FACTOR, COMPLETION_STRINGS, IN_FILE_PATHS, OUT_FILE_PATHS, SBATCH_POSTFIX
from settings import Settings
from slab import Slab, mol_bonded_to_slab
from molecule import Molecule
from filenames import *

import ase_custom

TEST = False #do not actually launch the jobs, simply prints the command


#OK (code agnostic) 
def is_completed(program : str, calc_type : str, completion_type : str, i_calc : int):
    '''
    Check if the given calculation is completed, reading the output file

    Args:
    - program: DFT program. Possible values: 'ESPRESSO' or 'VASP'
    - calc_type: 'SCREENING' or 'RELAX'
    - completion_type: 'RELAX_COMPLETED' or 'SCF_NONCONVERGED'
    - i_calc: numeric index of the calculation
    '''

    filename = OUT_FILE_PATHS[calc_type][program].format(i_calc)
    searchfor = COMPLETION_STRINGS[completion_type][program]
    
    with open(filename, 'r') as f:
        file_content = f.readlines()

    completed = False
    for line in file_content:
        if searchfor in line:
            completed = True
            break
    
    return completed

#OK (code agnostic) 
def _get_configurations_numbers():
    '''
    Function to get the list of indices from the site_labels.csv file

    Returns a list of indices
    '''

    site_labels = np.genfromtxt(labels_filename, delimiter=',', names=True)
    return site_labels['Label']

#OK (code agnostic) 
def _get_actually_present_outputs(program : str, calc_type : str):
    '''
    Function to get the list of indices of the outputs that are actually present.

    Returns a list of indices.
    '''
    
    indices = _get_configurations_numbers()

    existing_indices = []
    for index in indices:
        if os.path.isfile(OUT_FILE_PATHS[calc_type][program].format(index)):
            existing_indices.append(index)

    return existing_indices

#OK (code agnostic) 
def get_energy(program : str, calc_type : str, i_calc : int):
    '''
    Returns the TOTAL energy for a given configuration, or None if not available

    Args:
    - program: DFT program. Possible values: 'ESPRESSO' or 'VASP'
    - calc_type: 'SCREENING' or 'RELAX'
    - i_calc: numeric index of the calculation
    '''
    
    filename = OUT_FILE_PATHS[calc_type][program].format(i_calc)

    try:
        atoms = read(filename)
        return atoms.get_potential_energy()
    except:
        return None

#OK (code agnostic) 
def get_calculations_results(program : str, calc_type : str, E_slab_mol : list = [0,0], VERBOSE : bool =True):
    '''
    Returns a dictionary in the format

    {
        'energies': {1: -1200, 2: -1300},
        'relax_completed': {1: True, 2: False},
        'scf_nonconverged': {1: False, 2: False},
    }

    The calculations with no output file are not included in the dictionary

    Args:
    - program: DFT program. Possible values: 'ESPRESSO' or 'VASP'
    - calc_type: 'SCREENING' or 'RELAX'
    - E_slab_mol: if not [0,0] returns the adsorption energy in eV
    - VERBOSE: give warning messages for noncompleted calculations
    '''

    results = {'energies': {}, 'relax_completed': {}, 'scf_nonconverged': {}}

    indices = _get_configurations_numbers()

    for index in indices:
        if not os.path.isfile(OUT_FILE_PATHS[calc_type][program].format(index)):
            continue #skip if file does not exist yet

        energy = (get_energy(program, calc_type, index) - (E_slab_mol[0]+E_slab_mol[1])) * UNITS_TO_EV_FACTOR[program]
        relax_completed = is_completed(program, calc_type, 'RELAX_COMPLETED', index)
        scf_nonconverged = is_completed(program, calc_type, 'SCF_NONCONVERGED', index)

        results['energies'].update({index : energy})
        results['relax_completed'].update({index : relax_completed})
        results['scf_nonconverged'].update({index : scf_nonconverged})

        if(VERBOSE and scf_nonconverged): 
            print(f'Warning! Config. {index} failed to reach SCF convergence after electron_maxstep.\
                The last usable energy value will be used')

    
    return results

#OK (code agnostic) 
def write_results_to_file(TXT=False):
    '''
    Function to write the calculations results to a csv file.
    It can be called before all the jobs have finished.

    It needs to read the 'site_labels.csv' file.

    Args:
    - TXT: write a txt file (tab separated) instead of csv, sorted by screening or relax energies
    '''
    #

    settings = Settings()

    datafile = pd.read_csv(labels_filename, index_col=0)

    if os.path.isdir(screening_outdir): #for screening
        screening_results = \
            get_calculations_results(program=settings.program, calc_type='SCREENING', E_slab_mol=settings.E_slab_mol)
        
        column_name = 'Eads_scr(eV)' if 0 not in settings.E_slab_mol else 'Etot_scr(eV)'

        column_data = []
        for i in datafile.index: #if the file exists, and so the energy might be present, or it might be None
            if i in screening_results['energies']:
                
                column_data.append(screening_results['energies'][i] )

                if screening_results['scf_nonconverged'][i]:
                    print(f'Warning! {i} failed to reach SCF convergence after electron_maxstep.\
                           The energy will be marked with **')
                    column_data[-1] = str(column_data[-1])+'**'
                if not screening_results['relax_completed'][i]:
                    print(f'Warning! {i} relaxation has not reached final configuration. \
                           The energy will be marked with a *')
                    column_data[-1] = str(column_data[-1])+'*'

            else: #if the file does not exist
                column_data.append(None)

        datafile[column_name] = column_data


    if os.path.isdir(relax_outdir): #for relax
        relax_results = \
            get_calculations_results(program=settings.program, calc_type='RELAX', E_slab_mol=settings.E_slab_mol)
        
        column_name = 'Eads_rel(eV)' if 0 not in settings.E_slab_mol else 'Etot_rel(eV)'

        column_data = []
        for i in datafile.index: #if the file exists, and so the energy might be present, or it might be None
            if i in relax_results['energies']:
                
                column_data.append(relax_results['energies'][i] )

                if relax_results['scf_nonconverged'][i]:
                    print(f'Warning! {i} failed to reach SCF convergence after electron_maxstep.\
                           The energy will be marked with **')
                    column_data[-1] = str(column_data[-1])+'**'
                if not relax_results['relax_completed'][i]:
                    print(f'Warning! {i} relaxation has not reached final configuration. \
                           The energy will be marked with a *')
                    column_data[-1] = str(column_data[-1])+'*'

            else: #if the file does not exist
                column_data.append(None)

        datafile[column_name] = column_data

    
    #add bonding info
    slab = Slab(settings.slab_filename)
    mol = Molecule(settings.molecule_filename)
    mol_indices = np.arange(mol.natoms) if settings.mol_before_slab else np.arange(mol.natoms) + slab.natoms
    bonding_status = []
    for i in datafile.index: #if the file exists, and so the energy might be present, or it might be None
        if i in relax_results['energies']:
            status = check_bond_status(settings.program, 
                                                    calc_type= 'RELAX' if os.path.isdir(relax_outdir) else 'SCREENING', 
                                                    i_calc=i, 
                                                    mol_indices=mol_indices)
            bonding_status.append('Yes' if status else 'No')
        else: #if the file does not exist
            bonding_status.append(None)

    datafile['Bonded'] = column_data   
        

    if(TXT): #sort by energy column (relax if available, else screening)
        datafile.sort_values(by=column_name)
  
    datafile.to_csv(results_filename.replace('csv', 'txt' if TXT else 'csv'), sep='\t' if TXT else ',')

#OK (code agnostic) 
def check_bond_status(program : str, calc_type : str, i_calc : int, mol_indices : list):
    '''
    Reads output file and returns True if the molecule is bonded to the surface

    Args:
    - program: DFT program. Possible values: 'ESPRESSO' or 'VASP'
    - calc_type: 'SCREENING' or 'RELAX'
    - i_calc: numeric index of the calculation
    - mol_indices: indices of the atoms belonging to the molecule
    '''
    filename = OUT_FILE_PATHS[calc_type][program].format(i_calc)
    atoms = read(filename, results_required=False)

    slab = atoms.copy()
    del slab[[atom.index for atom in slab if atom.index in mol_indices]]
    mol = atoms.copy()
    del mol[[atom.index for atom in mol if atom.index not in mol_indices]]

    return mol_bonded_to_slab(slab, mol)

#OK (code agnostic) 
def launch_jobs(program : str, calc_type : str, jobscript : str, sbatch_command : str, indices_list : list):
    '''
    Launch the calculations.

    Args:
    - program: DFT program. Possible values: 'ESPRESSO' or 'VASP'
    - calc_type: 'SCREENING' or 'RELAX'
    - jobscript: path of the jobscript file
    - sbatch_command: command to submit the jobscript (in Slurm it is sbatch)
    - indices_list: indices of the calculations
    '''
    main_dir = os.getcwd()

    for index in indices_list:

        j_dir = f'{screening_outdir if calc_type is "SCREENING" else relax_outdir}/{index}'
        os.makedirs(j_dir, exist_ok=True) #for QE (for VASP they are already created by write_inputs)
        shutil.copyfile(jobscript, f'{j_dir}/{jobscript_stdname}')
        
        os.chdir(j_dir)   ####################
        
        #change job title (only for slumr jobscripts)
        with open(jobscript_stdname, 'rw') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if "job-name" in line:
                    lines[i] = f"{line.split('=')[0]} = {'scr' if calc_type is 'SCREENING' else 'RELAX'}_{index}\n"
                    break        
            f.writelines(lines)


        launch_string = f"{sbatch_command} {jobscript_stdname} {SBATCH_POSTFIX[calc_type][program].format(main_dir, i)}"
        if(TEST): print(launch_string)
        else: os.system(launch_string)  #launchs the jobscript in j_dir from j_dir
        os.chdir(main_dir) ####################

#OK (code agnostic) 
def restart_jobs(calc_type : str):
    '''
    Restart the uncompleted calculations

    Args:
    - calc_type: 'SCREENING' or 'RELAX'
    '''

    from settings import Settings
    settings = Settings()

    existing_indices = _get_actually_present_outputs(settings.program, calc_type)
    indices_to_restart = [index for index in existing_indices if is_completed(settings.program, calc_type, 'RELAX_COMPLETED', index)]

    #edit input files
    edit_files_for_restart(settings.program, calc_type, indices_to_restart)  

    #launch the calculations
    main_dir = os.getcwd()
    for index in indices_to_restart:
        j_dir = f'{screening_outdir if calc_type is "SCREENING" else relax_outdir}/{index}'
        os.chdir(j_dir)

        launch_string = f"{settings.sbatch_command} {jobscript_stdname} {SBATCH_POSTFIX[calc_type][settings.program].format(main_dir, index)}"
        if(TEST): print(launch_string)
        else: os.system(launch_string)  #launchs the jobscript in j_dir from j_dir 

        os.chdir(main_dir)

#OK (code agnostic)
def saveas(calc_type : str, i_or_f : str, saveas_format : str):
    '''
    Save all the configurations in a different format, e.g. xyz or cif.

    Args:
    - calc_type: 'screening' or 'relax'
    - i_or_f: initial or final coordinates of the relaxation (both for screening and full relax)
    - saveas_format: file format, e.g. xyz
    '''

    settings = Settings()

    if i_or_f is 'i':
        FILE_PATHS = IN_FILE_PATHS
    elif i_or_f == 'f':
        FILE_PATHS = OUT_FILE_PATHS
    else: 
        raise RuntimeError(f"Wrong arguments: passed '{calc_type} {i_or_f}', expected 'screening i/f' or 'relax i/f'")
    if calc_type != 'screening' and calc_type != 'relax':
        raise RuntimeError(f"Wrong argument: passed '{calc_type}', expected 'screening' or 'relax'")

    folder = f"{saveas_format}/{calc_type}"

    print(f"Saving files to {folder}...")
    os.makedirs(folder, exist_ok=True)
    
    indices = _get_configurations_numbers()
    for i in indices:
        if os.path.isfile(FILE_PATHS[calc_type.upper()][settings.program].format(i)):
            atoms = read(FILE_PATHS[calc_type.upper()][settings.program].format(i))
            if(saveas_format == 'xyz'):
                ase_custom.write_xyz_custom(f'{folder}/{calc_type}_{i}.{saveas_format}', atoms)
            else:
                write(f'{folder}/{calc_type}_{i}.{saveas_format}', atoms)

    print("All files saved.")
    