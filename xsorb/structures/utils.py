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

from xsorb.adsorptiondata.adsorptioncalculation import BondInfo


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


def slab_mol_bonds(slab : Atoms, mol: Atoms, mult : float = 1.1) -> list[BondInfo] | None:
    '''
    Returns a list of BondInfo objects representing the bonds between the slab and the molecule,
    or None if no bonds are found.
    Based on covalent radii (with mult factor of 1.1)

    Args:
    - slab: Atoms object for the slab
    - mol: Atoms object for the molecule
    - mult: factor to multiply the covalent radii

    Returns:
    - List of BondInfo objects if bonds are found, otherwise None
    '''

    atoms : Atoms = slab+mol
    cutoffs = natural_cutoffs(atoms, mult=mult)
    nl = NeighborList(cutoffs, skin=0, self_interaction=False, bothways=True)
    nl.update(atoms)
    cm = nl.get_connectivity_matrix()
    dm = atoms.get_all_distances(mic=True)

    bonds_list : list[BondInfo] = []
    for i in range(len(slab)):
        for j in range(len(mol)):
            if cm[i, len(slab)+j]:
                bonds_list.append(
                    BondInfo(
                        slab_atom_id=i,
                        mol_atom_id=j,
                        slab_atom_species=slab.get_chemical_symbols()[i],
                        mol_atom_species=mol.get_chemical_symbols()[j],
                        length=dm[i, len(slab)+j]
                    )
                )

    if len(bonds_list) == 0:
        return None
    else:
        return bonds_list
