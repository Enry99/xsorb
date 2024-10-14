#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue 28 Feb 2023
@author: Enrico Pedretti

Function definitions to read from pwo and launch scripts

"""

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
