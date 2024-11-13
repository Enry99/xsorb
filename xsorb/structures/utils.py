#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Enrico Pedretti

'''
Small module with utility functions for structures
'''

from __future__ import annotations

from ase import Atoms
from ase.constraints import FixCartesian
from ase.neighborlist import NeighborList, natural_cutoffs


def set_fixed_slab_constraints(atoms : Atoms, slab_indices : list | None = None) -> None:
    '''
    Inplace modifies the Atoms object to fully fix the slab atoms.
    '''
    indices = slab_indices if slab_indices is not None else list(range(len(atoms)))
    #get indices of already present constraints
    slab_constraints = [FixCartesian(atom_index) for atom_index in indices]
    mol_constraints = [constraint for constraint in atoms.constraints \
                       if constraint.index[0] not in slab_indices]
    atoms.set_constraint(slab_constraints + mol_constraints)


def slab_mol_bonds(slab : Atoms, mol: Atoms, mult : float = 1.1) -> str:
    '''
    Returns the list of the bonds between the molecule and the slab.
    Based on covalent radii (with mult factor of 1.15)

    Args:
    - slab: Atoms object for the slab
    - mol: Atoms object for the molecule
    - mult: factor to multiply the covalent radii

    Returns:
    - str, list of the bonds in the format 'Cu1-O2,Cu5-O3', or 'none' if no bonds are found
    '''

    atoms = slab+mol
    cutoffs = natural_cutoffs(atoms, mult=mult)
    nl = NeighborList(cutoffs, skin=0, self_interaction=False, bothways=True)
    nl.update(atoms)
    cm = nl.get_connectivity_matrix()

    bonds_list = []
    for i in range(len(slab)):
        for j in range(len(mol)):
            if cm[i, len(slab)+j]:
                bonds_list.append(f'{mol.symbols[j]}{j}-{slab.symbols[i]}{i}')

    if len(bonds_list) == 0:
        return 'none'
    else:
        return ','.join(bonds_list)
