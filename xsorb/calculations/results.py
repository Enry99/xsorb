#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Enrico Pedretti

'''
Module to check and update the results of the calculations
'''
from __future__ import annotations
from pathlib import Path
import sys
import logging

from xsorb.structures.utils import slab_mol_bonds
from xsorb.ase_custom.io import ase_custom_read as read
from xsorb.dft_codes.definitions import (
    SCF_NONCONVERGED_STRINGS, SCF_CONVERGED_STRINGS, OPTIMIZATION_COMPLETED_STRINGS)
from xsorb.settings.settings import Energies
from xsorb.adsorptiondata import AdsorptionCalculation
from xsorb.adsorptiondata.adsorptioncalculation import CalculationResults
from xsorb.io.utils import progressbar


def is_optimization_completed(filename : str, program : str):
    '''
    Check if the given calculation is completed, reading the output file

    Args:
    - filename: path to the LOG_FILE(==output file for espresso)
    - calc_type: 'screening','relax','mlopt'

    Returns:
    True or False
    '''

    searchfor = OPTIMIZATION_COMPLETED_STRINGS[Path(filename).suffix.lstrip('.')]

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
    - calc_type: 'screening','relax','mlopt'

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


def get_atoms_traj_from_calc(filename : str):
    '''
    Reads the output file and returns the atoms object.
    If there is an error reading the file, it prints a message and returns None

    Args:
    - filename: path to the output file

    Returns:
    list[Atoms] or None
    '''

    try:
        trajectory = read(filename, index=':')
        if not isinstance(trajectory, list):
            trajectory = [trajectory]
    except Exception as exc: #pylint: disable=broad-except
        logging.error(f'Error reading file {filename}: {exc}.')
        return None

    return trajectory


def get_bond_status(atoms, mol_indices : list, mult : float):
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


def update_calculations_results(*,systems: list[AdsorptionCalculation],
                             program : str,
                             energies : Energies,
                             mult : float,
                             store_full_trajectories : bool = True,
                             verbose : bool =True):
    '''
    Reads the output files and inplace updates the calculation results

    Args:
    - systems: list of AdsorptionCalculation objects
    - program: 'espresso','vasp','ml'
    - energies: Energies dataclass with slab and molecule energies
    - mult: multiplicative factor for the covalent radii to determine bonding.
    - store_full_trajectories: if True, saves the full trajectory in the CalculationResults
    - verbose: if True, prints warnings if files are missing
    '''

    for system in progressbar(systems, prefix='Reading calculation results: ',
                              transparent=not verbose):
        assert system.calc_info is not None #DEBUG

        if not Path(system.calc_info.out_file_path).exists() \
            or not Path(system.calc_info.log_file_path).exists():
            if verbose:
                if not Path(system.calc_info.out_file_path).exists():
                    missing_file = system.calc_info.out_file_path
                else:
                    missing_file = system.calc_info.log_file_path
                logging.warning(f'Warning! File {missing_file} not found. Skipping.')
            continue

        traj = get_atoms_traj_from_calc(system.calc_info.out_file_path)
        if not traj:
            continue

        try:
            # read results from file
            atoms = traj[-1]

            e_slab = energies.E_slab_ml if program == 'ml' else energies.E_slab
            e_mol = energies.E_mol_ml if program == 'ml' else energies.E_mol
            if isinstance(e_mol, list):
                e_mol = e_mol[system.adsorption_structure.mol_rot.conform_id]
            total_e_slab_mol = e_slab + e_mol

            adsorption_energy = atoms.get_potential_energy() - total_e_slab_mol

            adsorption_energy_evol = \
                    [at.get_potential_energy() - total_e_slab_mol for at in traj]

            if is_optimization_completed(system.calc_info.log_file_path, program):
                status = 'completed'
            elif is_scf_not_converged(system.calc_info.log_file_path, program):
                status = 'scf_nonconverged'
            else:
                status = 'incomplete'

            mol_indices = system.adsorption_structure.mol_indices
            bonds = get_bond_status(atoms, mol_indices, mult)

            mol_ref_idx = system.adsorption_structure.mol_rot.mol_atom
            adsize_z = system.adsorption_structure.adsite.coords[2]
            if mol_ref_idx == -1:
                #geometric center of the molecule
                final_dz = atoms.positions[mol_indices].mean(axis=0)[2] - adsize_z
            else:
                #reference atom of the molecule
                mol_ref_index = mol_ref_idx + mol_indices[0]
                final_dz = atoms[mol_ref_index].position[2] - adsize_z

            #inplace update of the calculation info and results
            system.calc_info.status = status

            system.calc_results = CalculationResults(
                atoms=atoms,
                adsorption_energy=adsorption_energy,
                bonds=bonds,
                trajectory=traj if store_full_trajectories else [atoms],
                adsorption_energy_evol=adsorption_energy_evol,
                final_dz=final_dz
            )


        except Exception as e: #pylint: disable=broad-except
            logging.info(f'Unable to read energy in  file {system.calc_info.out_file_path},'\
                        f' error: {e}. The calculation may be still running. Skipping.')
