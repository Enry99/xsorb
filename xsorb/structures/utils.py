import numpy as np
from ase import Atoms
from ase.geometry import get_distances
from ase.neighborlist import NeighborList, natural_cutoffs

def closest_pair(slab : Atoms, mol: Atoms):
    '''
    Returns the indices of the closest pair of (slab, molecule) atoms, and their distance

    Args:
    - slab: Atoms object for the slab
    - mol: Atoms object for the molecule 
    '''

    dist_matrix = get_distances(mol.positions, slab.positions, slab.cell, pbc=True)[1]

    i_mol, j_slab = np.unravel_index(np.argmin(dist_matrix), dist_matrix.shape)
    mindist = dist_matrix[i_mol, j_slab]

    return i_mol, j_slab, mindist


def slab_mol_bonds(slab : Atoms, mol: Atoms):
    '''
    Returns the list of the bonds between the molecule and the slab.

    Args:
    - slab: Atoms object for the slab
    - mol: Atoms object for the molecule
    '''

    atoms = slab+mol
    cutoffs = natural_cutoffs(atoms, mult=1.15)
    nl = NeighborList(cutoffs, skin=0, self_interaction=False, bothways=True)
    nl.update(atoms)
    cm = nl.get_connectivity_matrix()

    bonds_list = []
    for i in range(len(slab)):
        for j in range(len(mol)):
            if cm[i, len(slab)+j]:
                bonds_list.append('{0}{1}-{2}{3}'.format(mol.symbols[j], j, slab.symbols[i], i))
    
    return bonds_list
