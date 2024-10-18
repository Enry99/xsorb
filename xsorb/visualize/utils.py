'''
Utility function to center the molecule
'''

from pathlib import Path
import json

import numpy as np
from ase import Atoms

def read_custom_colors():

    if Path("custom_colors.json").exists():
        with open("custom_colors.json", "r") as f:
            custom_colors = json.load(f)
            print(f"Custom colors read from file.")
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


def plot_overview_grid(calc_type : str,
                       outfiles : list[str],
                       calc_indices : list[int],
                       energies : list[float],
                       stars : list[str]):
    '''
    Plot a grid with the images of the calculations,
    with the numbers (calc indices) and the energies.
    The lowest energy is highlighted in red.

    Args:
    - calc_type: 'initial','screening','relax','ml_opt'
    - outfiles: list of the paths to the images
    - calc_indices: list of the calculation indices
    - energies: list of the energies
    - stars: list of the stars ('*' if not completed, '**' if not converged)
    '''

    from matplotlib import pyplot as plt
    import matplotlib.image as mpimg

    #calculate aspect ratio for one image:
    img = mpimg.imread(outfiles[0])
    ar = float(img.shape[1]) / float(img.shape[0]) #width/height

    #calculate the number of rows and columns for the grid, and setup the figure
    n_rows_fig = max(int(np.ceil(len(calc_indices)/5.)), 1)
    n_cols_fig = max(int(np.ceil(len(calc_indices)/float(n_rows_fig))), 1)
    height = n_rows_fig * 2
    width  = n_cols_fig * ar * 2
    fig = plt.figure(figsize=(width, height))
    fig.subplots_adjust(wspace=0.1)
    axes = [fig.add_subplot(n_rows_fig,n_cols_fig,i) for i in range(1,len(calc_indices) + 1)]

    # find the index of the minimum energy
    imin = np.argmin(energies)
    # read the images from files, and plot them into the grid
    for i, calc_index in enumerate(calc_indices):

        img = mpimg.imread(outfiles[i])
        axes[i].imshow(img)
        axes[i].axis('equal') #ensures all images are in identical rectangular boxes
        axes[i].set_title(f'{energies[i]:.2f}{stars[i]} eV', fontsize = 5, pad=1,
                          color='red' if i == imin else 'black')

        #rect = plt.patches.Rectangle((0.0, 0.93), 0.08, 0.07, transform=axes[i].transAxes,
        # facecolor='white', edgecolor='black', linewidth=0.5)
        #axes[i].add_artist(rect)

        axes[i].text(0.045, 0.988, str(calc_index), fontsize = 3.5, transform=axes[i].transAxes,
                     horizontalalignment="left", verticalalignment="top",
                     bbox=dict(boxstyle='square', linewidth=0.5, fc="w", ec="k"))
        axes[i].set_xticks([])
        axes[i].set_yticks([])

    fig.savefig(f"{calc_type}_overview.png", dpi=700, bbox_inches='tight')
