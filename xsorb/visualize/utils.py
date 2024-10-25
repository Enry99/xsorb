#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
- Utility function to center the molecule
- Utility function to read custom colors from a file

'''

from __future__ import annotations
from pathlib import Path
import json
import sys

import numpy as np
from ase import Atoms

def read_custom_colors():
    '''
    Read the custom colors from the file custom_colors.json and return a dictionary.
    If the file does not exist, return None.
    '''

    if Path("custom_colors.json").exists():
        with open("custom_colors.json", "r", encoding=sys.getfilesystemencoding()) as f:
            custom_colors = json.load(f)
            print("Custom colors read from file.")
        return custom_colors

    return None


def get_centered_mol_and_slab(atoms : Atoms, mol_indices : list[int]):
    """
    Get the centered molecule and slab.

    Parameters:
    - atoms (Atoms): The atoms object with the molecule and the slab
    - mol_indices (list): The indices of the atoms in the molecule

    Returns:
    atoms (Atoms): The atoms with centered molecule
    mol_indices (list): The indices of the atoms in the molecule
    transl_vector (np.array): The translation vector to center the molecule
    """

    from ase import neighborlist
    from scipy import sparse

    #avoid modifying the original structure
    atoms = atoms.copy()
    mol = atoms[mol_indices]
    slab = atoms[[i for i in range(len(atoms)) if i not in mol_indices]]

    #replicate in both directions to have at least one fully connected molecule
    mol.set_constraint() #suppress warnings
    mol_replicated = mol * [2,2,1]
    mol_replicated.set_pbc(False)

    #get the indices of the atoms in each molecule (connected components)
    cutoffs = neighborlist.natural_cutoffs(mol_replicated, mult=1.2)
    nl = neighborlist.NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(mol_replicated)
    cmatrix = nl.get_connectivity_matrix()
    n_components, component_list = sparse.csgraph.connected_components(cmatrix)
    molecules_indices = [ [ atom_id for atom_id in range(len(component_list)) \
                           if component_list[atom_id] == molecule_id ] \
                            for molecule_id in range(n_components) ]

    #find the molecule with the highest number of atoms (which is therefore fully connected)
    max_molecule_idx = np.argmax([len(molecule) for molecule in molecules_indices])
    mol_replicated = mol_replicated[molecules_indices[max_molecule_idx]]

    #fin the geometric center of the fully connected molecule,
    # that will be translated at the center of the cell
    mol_center = mol_replicated.positions.mean(axis=0)
    mol_center[2] = 0 #we don't want to translate in the z direction (just xy)

    cell_center = slab.cell[:][0]/2 + slab.cell[:][1]/2  #a1/2 + a2/2 (vector sum)

    transl_vector = cell_center - mol_center

    mol.translate(transl_vector)
    slab.translate(transl_vector)
    #wrap both since we translated them (also the molecule, since
    # we might have chosen a different replica)
    mol.wrap()
    slab.wrap()

    atoms = slab + mol
    new_mol_indices = list(range(len(slab), len(atoms)))

    return atoms, new_mol_indices, transl_vector
