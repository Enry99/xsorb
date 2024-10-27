#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue 28 Feb 2023
@author: Enrico Pedretti

Function definitions to read from pwo and launch scripts

"""

from __future__ import annotations
import os, shutil, subprocess
import numpy as np
import pandas as pd
from ase.io import read, write
from ase.calculators.calculator import PropertyNotImplementedError
from xsorbed.dftcode_specific import edit_files_for_restart, OPTIMIZATION_COMPLETED_STRINGS, SCF_NONCONVERGED_STRINGS, SCF_CONVERGED_STRINGS, IN_FILE_PATHS, OUT_FILE_PATHS, LOG_FILE_PATHS, SBATCH_POSTFIX
from xsorbed.settings import Settings
from xsorbed.slab import Slab, slab_mol_bonds
from xsorbed.molecule import Molecule
from xsorbed.common_definitions import *

from xsorbed import ase_custom

TEST = False #do not actually launch the jobs, simply prints the command



def optimization_completed(program : str, calc_type : str, i_calc : int):
    '''
    Check if the given calculation is completed, reading the output file

    Args:
    - program: DFT program. Possible values: 'ESPRESSO' or 'VASP'
    - calc_type: 'SCREENING' or 'RELAX'
    - i_calc: numeric index of the calculation

    Returns:
    True or False
    '''

    filename = LOG_FILE_PATHS[calc_type][program].format(i_calc)
    searchfor = OPTIMIZATION_COMPLETED_STRINGS[program]

    with open(filename, 'r') as f:
        file_content = f.readlines()

    completed = False
    for line in file_content:
        if searchfor in line:
            completed = True
            break

    return completed

def scf_not_converged(program : str, calc_type : str, i_calc : int):
    '''
    Check if the given calculation has not reached SCF convergence, reading the output file

    Args:
    - program: DFT program. Possible values: 'ESPRESSO' or 'VASP'
    - calc_type: 'SCREENING' or 'RELAX'
    - i_calc: numeric index of the calculation

    Returns:
    True or False
    '''

    filename = OUT_FILE_PATHS[calc_type][program].format(i_calc)
    searchfor = SCF_NONCONVERGED_STRINGS[program]
    convergence_string = SCF_CONVERGED_STRINGS[program]

    with open(filename, 'r') as f:
        file_content = f.readlines()

    # the last one (conv or not conv) determines the status
    nonconv = False
    for line in file_content:
        if searchfor in line:
            nonconv = True
            continue
        if convergence_string in line:
            nonconv = False

    return nonconv


def _get_configurations_numbers():
    '''
    Function to get the list of indices from the site_labels.csv file

    Returns a list of indices
    '''

    site_labels = np.genfromtxt(labels_filename, delimiter=',', names=True, dtype=int)
    return site_labels['Label']


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


def get_energy(program : str, calc_type : str, i_calc : int, full_evolution : bool = False, VERBOSE : bool = True):
    '''
    Returns the TOTAL energy for a given configuration, or None if not available

    Args:
    - program: DFT program. Possible values: 'ESPRESSO' or 'VASP'
    - calc_type: 'SCREENING' or 'RELAX'
    - i_calc: numeric index of the calculation
    - full_evolution: if True, returns a NDarray with the energies of the relaxation
    '''

    filename = OUT_FILE_PATHS[calc_type][program].format(i_calc)

    traj = read(filename, index=':')
    energies = []

    for atoms in traj:
        try:
            energy = atoms.get_potential_energy()
        except PropertyNotImplementedError as e:
            #if VERBOSE: print(f'Warning: unable to read last energy value from file {filename}.')
            continue
        energies.append(energy)

    if len(energies) == 0:
        if VERBOSE: print(f'Warning: unable to read ANY energy from file {filename}.')
        return None
    elif full_evolution:
        return np.array(energies)
    else:
        return energies[-1]

def get_energy_ml(filename : str, full_evolution : bool = False):

    try:
        if not full_evolution:
            atoms = read(filename)
            return atoms.get_potential_energy()
        else:
            atoms_list = read(filename, index=':')
            return np.array([a.get_potential_energy() for a in atoms_list])
    except:
        print(f'Warning: unable to read energy from file {filename}.')
        return None


def get_calculations_results(program : str,
                             calc_type : str,
                             E_slab_mol : list = [0,0],
                             full_evolution : bool = False,
                             VERBOSE : bool =True):
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
    - full_evolution: if True, each config in 'energies' will contain an array with the energies during relaxation
    - VERBOSE: give warning messages for noncompleted calculations
    '''

    results = {'energies': {}, 'relax_completed': {}, 'scf_nonconverged': {}}

    indices = _get_configurations_numbers()

    for index in indices:
        if not os.path.isfile(OUT_FILE_PATHS[calc_type][program].format(index)):
            continue #skip if file does not exist yet

        energy = get_energy(program, calc_type, index, full_evolution, VERBOSE=VERBOSE)


        if energy is not None:
            energy -= (E_slab_mol[0]+E_slab_mol[1])
        relax_completed = optimization_completed(program, calc_type, index)
        scf_nonconverged = scf_not_converged(program, calc_type, index)

        results['energies'].update({index : energy})
        results['relax_completed'].update({index : relax_completed})
        results['scf_nonconverged'].update({index : scf_nonconverged})

        if(VERBOSE and scf_nonconverged):
            print(f'Warning! Config. {index} failed to reach SCF convergence after electron_maxstep. The last usable energy value will be used')


    return results


def get_calculations_results_ml(full_evolution : bool = False, subtract_eslabmol : bool = True, return_separately_eslabmol : bool = False):
    '''
    Returns a dictionary in the format

    {
        'energies': {1: -1200, 2: -1300},
        'relax_completed': {1: True, 2: False},
    }

    The calculations with no output file are not included in the dictionary

    Args:
    - full_evolution: if True, each config in 'energies' will contain an array with the energies during relaxation
    - VERBOSE: give warning messages for noncompleted calculations
    '''

    results = {'energies': {}, 'relax_completed': {}}

    indices = _get_configurations_numbers()


    # extract e_slab and e_mol

    if subtract_eslabmol or return_separately_eslabmol:
        try:
            eslab = get_energy_ml(preopt_outdir+'/slab/slab.traj')
            emol = get_energy_ml(preopt_outdir+'/mol/mol.traj')
            E_slab_mol = [eslab, emol]
        except:
            eslab = None
            emol = None
            print('Warning: unable to read slab/molecule energies from preopt files. Total energies will be displayed.')

    for index in indices:
        if not os.path.isfile(OUT_FILE_PATHS['PREOPT']['ML'].format(index)):
            continue #skip if file does not exist yet

        energy = get_energy_ml(OUT_FILE_PATHS['PREOPT']['ML'].format(index), full_evolution)

        if energy is not None and subtract_eslabmol and E_slab_mol[0]:
            energy -= (E_slab_mol[0]+E_slab_mol[1])
        relax_completed = optimization_completed('ML', 'PREOPT', index)

        results['energies'].update({index : energy})
        results['relax_completed'].update({index : relax_completed})

    if return_separately_eslabmol:
        return results, eslab, emol
    else:
        return results



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


    if os.path.isdir(preopt_outdir): #for preopt
            preopt_results = \
                get_calculations_results_ml()

            column_name = 'Eads_pre(eV)'

            column_data = []
            for i in datafile.index: #if the file exists, and so the energy might be present, or it might be None
                if i in preopt_results['energies']:
                    if preopt_results['energies'][i] is None:
                        print(f'Config. {i} has no energy. It will be skipped.')
                        column_data.append(None)
                        continue

                    column_data.append(preopt_results['energies'][i] )

                    if not preopt_results['relax_completed'][i]:
                        print(f'Warning! {i} relaxation has not reached final configuration. The energy will be marked with a *')
                        column_data[-1] = f'{column_data[-1]:.3f}*'

                else: #if the file does not exist
                    column_data.append(None)

            datafile[column_name] = column_data


    if os.path.isdir(screening_outdir): #for screening
        screening_results = \
            get_calculations_results(program=settings.program, calc_type='SCREENING', E_slab_mol=settings.E_slab_mol, VERBOSE=False)

        column_name = 'Eads_scr(eV)' if 0 not in settings.E_slab_mol else 'Etot_scr(eV)'

        column_data = []
        for i in datafile.index: #if the file exists, and so the energy might be present, or it might be None
            if i in screening_results['energies']:
                if screening_results['energies'][i] is None:
                    print(f'Config. {i} has not reached the first scf convergence. It will be skipped.')
                    column_data.append(None)
                    continue

                column_data.append(screening_results['energies'][i] )

                if screening_results['scf_nonconverged'][i]:
                    print(f'Warning! {i} failed to reach SCF convergence in the last ionic step. The energy will be marked with **')
                    column_data[-1] = f'{column_data[-1]:.3f}**'
                elif not screening_results['relax_completed'][i]:
                    print(f'Warning! {i} relaxation has not reached final configuration. The energy will be marked with a *')
                    column_data[-1] = f'{column_data[-1]:.3f}*'

            else: #if the file does not exist
                column_data.append(None)

        datafile[column_name] = column_data


    if os.path.isdir(relax_outdir): #for relax
        relax_results = \
            get_calculations_results(program=settings.program, calc_type='RELAX', E_slab_mol=settings.E_slab_mol, VERBOSE=False)

        column_name = 'Eads_rel(eV)' if 0 not in settings.E_slab_mol else 'Etot_rel(eV)'

        column_data = []
        for i in datafile.index: #if the file exists, and so the energy might be present, or it might be None
            if i in relax_results['energies']:
                if relax_results['energies'][i] is None:
                    print(f'Config. {i} has not reached the first scf convergence. It will be skipped.')
                    column_data.append(None)
                    continue

                column_data.append(relax_results['energies'][i] )

                if relax_results['scf_nonconverged'][i]:
                    print(f'Warning! {i} failed to reach SCF convergence in the last ionic step. The energy will be marked with **')
                    column_data[-1] = f'{column_data[-1]:.3f}**'
                elif not relax_results['relax_completed'][i]:
                    print(f'Warning! {i} relaxation has not reached final configuration. The energy will be marked with a *')
                    column_data[-1] = f'{column_data[-1]:.3f}*'

            else: #if the file does not exist
                column_data.append(None)

        datafile[column_name] = column_data

    #add bonding info
    slab = Slab(settings.slab_filename)
    mol = Molecule(settings.molecule_filename, atom_index=settings.selected_atom_index)
    mol_indices = np.arange(mol.natoms) if settings.mol_before_slab else np.arange(mol.natoms) + slab.natoms
    bonding_status = []
    for i in datafile.index: #if the file exists, and so the energy might be present, or it might be None
        if os.path.isdir(relax_outdir) and i in relax_results['energies'] and relax_results['energies'][i]: #prioritize status from relax over screening
            status = check_bond_status(settings.program,
                                                    calc_type='RELAX',
                                                    i_calc=i,
                                                    mol_indices=mol_indices)
            bonding_status.append(','.join(status) if status else 'No')
        elif os.path.isdir(screening_outdir) and i in screening_results['energies'] and screening_results['energies'][i]:
            status = check_bond_status(settings.program,
                                                    calc_type='SCREENING',
                                                    i_calc=i,
                                                    mol_indices=mol_indices)
            bonding_status.append(','.join(status) if status else 'No')
        elif os.path.isdir(preopt_outdir) and i in preopt_results['energies'] and preopt_results['energies'][i]:
            status = check_bond_status('ML',
                                                    calc_type='PREOPT',
                                                    i_calc=i,
                                                    mol_indices=mol_indices)
            bonding_status.append(','.join(status) if status else 'No')
        else: #if the file does not exist
            bonding_status.append(None)

    datafile['Bonded'] = bonding_status


    if(TXT): #sort by energy column (relax if available, else screening)
        datafile.sort_values(by=column_name)

    datafile.to_csv(results_filename.replace('csv', 'txt' if TXT else 'csv'), sep='\t' if TXT else ',')

    print('Results file written.')


def check_bond_status(program : str, calc_type : str, i_calc : int, mol_indices : list):
    '''
    Reads output file and returns the list of bonds between slab and molecule

    Args:
    - program: DFT program. Possible values: 'ESPRESSO' or 'VASP'
    - calc_type: 'SCREENING' or 'RELAX'
    - i_calc: numeric index of the calculation
    - mol_indices: indices of the atoms belonging to the molecule
    '''
    filename = OUT_FILE_PATHS[calc_type][program].format(i_calc)
    atoms = read(filename)

    slab = atoms[[atom.index for atom in atoms if atom.index not in mol_indices]]
    mol = atoms[mol_indices]

    return slab_mol_bonds(slab, mol)


def launch_jobs(program : str, calc_type : str, jobscript : str, sbatch_command : str, indices_list : list, jobname_prefix : str = ''):
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

    submitted_jobs = []

    for index in indices_list:

        j_dir = f'{screening_outdir if calc_type == "SCREENING" else relax_outdir}/{index}'
        os.makedirs(j_dir, exist_ok=True) #for QE (for VASP they are already created by write_inputs)
        shutil.copyfile(jobscript, f'{j_dir}/{jobscript_stdname}')

        os.chdir(j_dir)   ####################

        #change job title (only for slumr jobscripts)
        with open(jobscript_stdname, 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if "job-name" in line:
                    prefix = jobname_prefix + ('_' if len(jobname_prefix) else '')
                    lines[i] = f"{line.split('=')[0]}={prefix}{'s' if calc_type == 'SCREENING' else 'r'}{index}\n"
                    break
        with open(jobscript_stdname, 'w') as f:
            f.writelines(lines)


        launch_string = f"{sbatch_command} {jobscript_stdname} {SBATCH_POSTFIX[calc_type][program].format(main_dir, index)}"
        if(TEST): print(launch_string)
        else:
            outstring = subprocess.getoutput(launch_string)  #launches the jobscript in j_dir from j_dir
            print(outstring)
            submitted_jobs.append(outstring.split()[-1])
        os.chdir(main_dir) ####################

    with open("submitted_jobs.txt", "a") as f:
        f.write("\n".join(submitted_jobs)+'\n')


def launch_jobs_ml(jobscript_ml : str,
                   sbatch_command : str,
                   indices_list : list = None,
                   explicit_labels : list | str = None,
                   fix_bondlengths : bool = False,
                   fix_slab : bool = False,
                   slab_indices : list = None,
                   jobname_prefix : str = ''):
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

    import pathlib
    xsorb_dir = pathlib.Path(__file__).parent.resolve()

    submitted_jobs = []

    if explicit_labels:
        if type(explicit_labels) == str:
            explicit_labels = [explicit_labels]
        for label in explicit_labels:
            j_dir = preopt_outdir+'/{0}'.format(label)
            shutil.copyfile(jobscript_ml, f'{j_dir}/{jobscript_stdname}')

            #change job title (only for slumr jobscripts)
            with open(f'{j_dir}/{jobscript_stdname}', 'r') as f:
                lines = f.readlines()
                for i, line in enumerate(lines):
                    if "job-name" in line:
                        lines[i] = f"{line.split('=')[0]}={jobname_prefix}_{'pre'}{label}\n"
                        break
            with open(f'{j_dir}/{jobscript_stdname}', 'w') as f:
                f.writelines(lines)

            os.chdir(j_dir)
            in_file = label+'.xyz'
            out_file = label+'.traj'
            log_file = label+'.log'

            if fix_bondlengths and fix_slab:
                fixstring = 'fixslab+fixmolbonds'
            elif fix_bondlengths:
                fixstring = 'fixbonds'
            elif fix_slab and label != 'mol':
                fixstring = 'fixslab'
            else:
                fixstring = 'nofix'

            launch_string = f"{sbatch_command} {jobscript_stdname} {xsorb_dir}/ml_opt.py {in_file} {out_file} {log_file} {main_dir} \
                {fixstring} {slab_indices[0]} {slab_indices[1]}"
            if(TEST): print(launch_string)
            else:
                outstring = subprocess.getoutput(launch_string)
                print(outstring)
                submitted_jobs.append(outstring.split()[-1])
            os.chdir(main_dir)
    else:

        for index in indices_list:

            j_dir = preopt_outdir+'/{0}'.format(index)
            shutil.copyfile(jobscript_ml, f'{j_dir}/{jobscript_stdname}')

            #change job title (only for slumr jobscripts)
            with open(f'{j_dir}/{jobscript_stdname}', 'r') as f:
                lines = f.readlines()
                for i, line in enumerate(lines):
                    if "job-name" in line:
                        prefix = jobname_prefix + ('_' if len(jobname_prefix) else '')
                        lines[i] = f"{line.split('=')[0]}={prefix}{'p'}{index}\n"
                        break
            with open(f'{j_dir}/{jobscript_stdname}', 'w') as f:
                f.writelines(lines)

            os.chdir(j_dir)
            in_file = os.path.basename(IN_FILE_PATHS['PREOPT']['ML'].format(index))
            out_file = os.path.basename(OUT_FILE_PATHS['PREOPT']['ML'].format(index))
            log_file = os.path.basename(LOG_FILE_PATHS['PREOPT']['ML'].format(index))

            if fix_bondlengths and fix_slab:
                fixstring = 'fixslab+fixmolbonds'
            elif fix_bondlengths:
                fixstring = 'fixbonds'
            elif fix_slab:
                fixstring = 'fixslab'
            else:
                fixstring = 'nofix'

            launch_string = f"{sbatch_command} {jobscript_stdname} {xsorb_dir}/ml_opt.py {in_file} {out_file} {log_file} {main_dir} \
                {fixstring} {slab_indices[0]} {slab_indices[1]}"
            if(TEST): print(launch_string)
            else:
                outstring = subprocess.getoutput(launch_string)
                print(outstring)
                submitted_jobs.append(outstring.split()[-1])
            os.chdir(main_dir)

    with open("submitted_jobs.txt", "a") as f:
        f.write("\n".join(submitted_jobs)+'\n')




def restart_jobs(calc_type : str):
    '''
    Restart the uncompleted calculations

    Args:
    - calc_type: 'SCREENING' or 'RELAX'
    '''

    from xsorbed.settings import Settings
    settings = Settings(read_energies=False)

    existing_indices = _get_actually_present_outputs(settings.program, calc_type)
    indices_to_restart = [index for index in existing_indices if not optimization_completed(settings.program, calc_type, index)]

    #edit input files
    edit_files_for_restart(settings.program, calc_type, indices_to_restart)

    #launch the calculations
    main_dir = os.getcwd()
    submitted_jobs = []
    for index in indices_to_restart:
        j_dir = f'{screening_outdir if calc_type == "SCREENING" else relax_outdir}/{index}'
        os.chdir(j_dir)

        launch_string = f"{settings.sbatch_command} {jobscript_stdname} {SBATCH_POSTFIX[calc_type][settings.program].format(main_dir, index)}"
        if(TEST): print(launch_string)
        else:
            outstring = subprocess.getoutput(launch_string)
            print(outstring)
            submitted_jobs.append(outstring.split()[-1])

        os.chdir(main_dir)

    with open("submitted_jobs.txt", "a") as f:
        f.write("\n".join(submitted_jobs)+'\n')


def scancel_jobs():
    '''
    Cancel all the jobs in the submitted_jobs.txt file, if they are running. Works only for slurm.
    '''

    if not os.path.isfile("submitted_jobs.txt"):
        print("submitted_jobs.txt file not found. No jobs to cancel.")
        return

    with open("submitted_jobs.txt", "r") as f:
        submitted_jobs = f.readlines()
        submitted_job_ids = [job.strip() for job in submitted_jobs]

    running_jobs = subprocess.getoutput("squeue --me").split("\n")[1:]
    running_job_ids = [job.split()[0] for job in running_jobs]

    for job in submitted_job_ids:
        if job in running_job_ids:
            print(f"Cancelling job {job}")
            os.system(f"scancel {job}")

    print("All jobs cancelled.")


def saveas(calc_type : str, i_or_f : str, saveas_format : str):
    '''
    Save all the configurations in a different format, e.g. xyz or cif.

    Args:
    - calc_type: 'screening' or 'relax'
    - i_or_f: initial or final coordinates of the relaxation (both for screening and full relax)
    - saveas_format: file format, e.g. xyz
    '''

    settings = Settings(read_energies=False)

    if i_or_f == 'i':
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
