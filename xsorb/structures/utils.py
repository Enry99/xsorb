'''
Small module with utility functions for structures
'''


from ase import Atoms
from ase.neighborlist import NeighborList, natural_cutoffs


def slab_mol_bonds(slab : Atoms, mol: Atoms):
    '''
    Returns the list of the bonds between the molecule and the slab.
    Based on covalent radii (with mult factor of 1.15)

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
                bonds_list.append(f'{mol.symbols[j]}{j}-{slab.symbols[i]}{i}')

    return bonds_list
