#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Enrico Pedretti

'''
Module to check and update the results of the calculations
'''
from __future__ import annotations
from typing import TYPE_CHECKING
from dataclasses import dataclass
from pathlib import Path
import sys

from ase import Atoms

import xsorb.io.jobs
from xsorb.structures.utils import slab_mol_bonds
from xsorb.ase_custom.io import ase_custom_read as read
from xsorb.dft_codes.definitions import (OUT_FILE_PATHS,
    SCF_NONCONVERGED_STRINGS, SCF_CONVERGED_STRINGS, OPTIMIZATION_COMPLETED_STRINGS)
if TYPE_CHECKING:
    from xsorb.io.inputs import WrittenSystem

@dataclass
class CalculationResults:
    '''
    Dataclass to store the results of a calculation
    '''
    atoms: Atoms
    adsorption_energy: float
    status : str #'completed', 'incomplete'
    scf_nonconverged : bool
    bonds : str
    job_status : str
    trajectory : list[Atoms]
    adsorption_energy_evol: list[float]
    final_dz: float


def is_optimization_completed(filename : str, program : str):
    '''
    Check if the given calculation is completed, reading the output file

    Args:
    - filename: path to the LOG_FILE(==output file for espresso)
    - calc_type: 'screening','relax','ml_opt'

    Returns:
    True or False
    '''

    searchfor = OPTIMIZATION_COMPLETED_STRINGS[program]

    with open(filename, 'r',encoding=sys.getfilesystemencoding()) as f:
        file_content = f.readlines()

    completed = False
    for line in file_content:
        if searchfor in line:
            completed = True
            break

    return completed


def is_scf_not_converged(filename : str, program : str):
    '''
    Check if the given calculation is completed, reading the output file

    Args:
    - filename: path to the LOG_FILE(==output file for espresso)
    - calc_type: 'screening','relax','ml_opt'

    Returns:
    True or False
    '''

    if program == 'ml': return False #pylint: disable=multiple-statements

    searchfor = SCF_NONCONVERGED_STRINGS[program]
    convergence_string = SCF_CONVERGED_STRINGS[program]

    with open(filename, 'r',encoding=sys.getfilesystemencoding()) as f:
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


def get_atoms_from_calc(filename : str, return_trajectory : bool = True):
    '''
    Reads the output file and returns the atoms object.
    If there is an error reading the file, it prints a message and returns None

    Args:
    - filename: path to the output file
    - return_trajectory: if True, returns a list of Atoms objects

    Returns:
    The atom object from the output file, or a list of Atoms objects if return_trajectory is True
    '''

    if return_trajectory:
        try:
            trajectory = read(filename, index=':')
        except Exception as exc:
            print(f'Error reading trajectory from file {filename}: {exc}. '\
                'Attempting to read only the last configuration.')
            try:
                trajectory = [read(filename)]
            except Exception as exc2:
                print(f'Error reading file {filename}: {exc2}.')
                return None

        return trajectory

    else:
        try:
            atoms = read(filename)
        except Exception as exc:
            print(f'Error reading file {filename}: {exc}.')
            return None

        return atoms


def bond_status(atoms, mol_indices : list, mult : float):
    '''
    Check the bonding status between the slab and the molecule.

    Args:
    - atoms: Atoms object
    - mol_indices: indices of molecule atoms in the slab+mol Atoms object
    - mult: multiplicative factor for the covalent radii to determine bonding.
    '''

    slab = atoms[[atom.index for atom in atoms if atom.index not in mol_indices]]
    mol = atoms[mol_indices]

    return slab_mol_bonds(slab, mol, mult)


def get_calculations_results(*,systems: list[WrittenSystem],
                             program : str,
                             total_e_slab_mol : float | None,
                             mult : float,
                             verbose : bool =True):
    '''
    Reads the output files and returns a list of CalculationResults objects.
    The calculations with no output will be added as None

    Args:
    - systems: list of WrittenSystem objects
    - program: DFT program. Possible values: 'espresso','vasp','ml'
    - total_e_slab_mol: total energy of the slab and molecule, if available
    - mult: multiplicative factor for the covalent radii to determine bonding.
    '''

    running_jobs = xsorb.io.jobs.get_running_jobs()

    results_list : list[CalculationResults | None ] = []

    for system in systems:
        if not Path(system.out_file_path).exists() or not Path(system.log_file_path).exists():
            if verbose:
                missing_file = system.out_file_path \
                    if not Path(system.out_file_path).exists() else system.log_file_path
                print(f'Warning! File {missing_file} not found. Skipping.')
            results_list.append(None)
            continue

        traj = get_atoms_from_calc(system.out_file_path)
        if not traj:
            results_list.append(None)
            continue

        try:
            atoms = traj[-1]
            adsorption_energy = atoms.get_potential_energy() - total_e_slab_mol

            adsorption_energy_evol = \
                    [at.get_potential_energy() - total_e_slab_mol for at in traj]

            status = 'completed' if is_optimization_completed(system.log_file_path, program) \
                else 'incomplete'

            scf_nonconverged = is_scf_not_converged(system.log_file_path, program)

            mol_indices = system.adsorption_structure.mol_indices
            bonds = bond_status(atoms, mol_indices, mult)

            mol_ref_idx = system.adsorption_structure.mol_rot.mol_atom
            adsize_z = system.adsorption_structure.adsite.coords[2]
            if mol_ref_idx == -1:
                #geometric center of the molecule
                final_dz = atoms.positions[mol_indices].mean(axis=0) - adsize_z
            else:
                #reference atom of the molecule
                mol_ref_index = mol_ref_idx + mol_indices[0]
                final_dz = atoms[mol_ref_index].position[2] - adsize_z

            job_status = 'running' if system.job_id in running_jobs else 'terminated'

            calc_results = CalculationResults(atoms=atoms,
                                            adsorption_energy=adsorption_energy,
                                            status=status,
                                            scf_nonconverged=scf_nonconverged,
                                            bonds=bonds,
                                            job_status=job_status,
                                            trajectory=traj,
                                            adsorption_energy_evol=adsorption_energy_evol,
                                            final_dz=final_dz)

            results_list.append(calc_results)

        except Exception as e:
            print(f'No energy in file {system.out_file_path}: ' \
                    f'possibly the calculation is still running. Error:{e}. Skipping.')
            results_list.append(None)

    return results_list
