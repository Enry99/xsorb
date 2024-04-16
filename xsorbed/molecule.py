#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue 28 Feb 2023
@author: Enrico Pedretti

Small helper class to manage the molecule

"""

import numpy as np
from matplotlib import pyplot as plt
from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core.bonds import CovalentBond
from ase.visualize.plot import plot_atoms
from ase.constraints import FixCartesian

from xsorbed import ase_custom 

#NOTE: ALL MODIFICATIONS DONE
class Molecule:
    '''
        Class to read molecule from file (e.g. Quantum ESPRESSO pwi/pwo or VASP POSCAR),
        find adsorption sites on the surface, and generate the adsorption structures with a molecule
        by placing the molecule on all the different sites

        Initialization parameters:
        - molecule_filename: file containing the structure of the molecule
        - atom_index: index of the reference atom (indexing starting from 0, in the order of the input file)
        - molecule_axis_atoms: indices of the two atoms defining the x-axis of the molecule (a = r2 - r1) 
        - axis_vector: [ax, ay, az] of the vector that is considered as the x-axis of the molecule
        - atoms_subset: indices of the atoms to include. Only removes atom NOT in this list, does not check
        if the atoms in the list actually exist in the molecule
        - break_bond_indices: indices of the two atoms of the bond to be broken. The fragment containing the first atom will be kept.
        - fixed_indices_mol: list of specific atoms to be fixed 
        (indices start from 0, with the ordering of the atoms in the input file)
        - fix_mol_xyz: which coordinates to fix for the fixed atoms, e.g. [True, True, False] = fix motion in x,y, free to move along z.
    '''


    def __init__(self, molecule_filename : str,
                 atom_index: int,
                 molecule_axis_atoms : list = None, 
                 axis_vector: list = None, 
                 atoms_subset : list = None,
                 break_bond_indices : list = None,
                 fixed_indices_mol : list = None, 
                 fix_mol_xyz : list = None):


        self.mol_ase = ase_custom.Atoms_custom(read(molecule_filename))

        #align axis to x axis
        if molecule_axis_atoms and axis_vector: raise RuntimeError("molecule axis cannot be given simultaneously as vector and by two atoms.")
        if molecule_axis_atoms:
            x1 = np.array(self.mol_ase.get_positions()[molecule_axis_atoms[0]])
            x2 = np.array(self.mol_ase.get_positions()[molecule_axis_atoms[1]])
            self.mol_ase.rotate(x2-x1, 'x')
        elif axis_vector:
            self.mol_ase.rotate(axis_vector, 'x')

        #SET CONSTRAINTS
        if fixed_indices_mol:
            c = [FixCartesian(idx, mask=[not x for x in fix_mol_xyz]) for idx in fixed_indices_mol]  #we need to negate: in qe 0 = fix, here 1(true)=fix
            self.mol_ase.set_constraint(c)
        else: self.mol_ase.set_constraint() #clean possible constraints read from file
        ###############################################################
            
        #Translate reference atom to origin
        if atom_index != -1:
            self.mol_ase.translate(-self.mol_ase.get_positions()[atom_index])
        else:
            self.mol_ase.translate(-self.mol_ase.get_center_of_mass())
            
        #select atoms subset
        if atoms_subset:            
            self.mol_ase = self.mol_ase[atoms_subset]

        elif break_bond_indices:
            molcopy = self.mol_ase.copy() 
            del molcopy.constraints #to suppress the warning about constraints not supported in pymatgen
            mol_pymat = AseAtomsAdaptor.get_molecule(molcopy)
            
            if not CovalentBond.is_bonded(mol_pymat[break_bond_indices[0]], mol_pymat[break_bond_indices[1]]):
                raise ValueError('The two selected atoms to split the molecule are not bonded.')
            
            mol_pymat = mol_pymat.break_bond(break_bond_indices[0], break_bond_indices[1])[0]

            included_indices = []
            for mol_atom in self.mol_ase:
                for frag_atom in mol_pymat:
                    if np.allclose(mol_atom.position, frag_atom.coords):
                        included_indices.append(mol_atom.index)
                        break

            self.mol_ase = self.mol_ase[included_indices]

        if atom_index != -1:
            ref_idx = [atom.index for atom in self.mol_ase if np.allclose(atom.position, [0,0,0])]
            if len(ref_idx) != 1: raise ValueError('Reference atom index not found.')
            self.reference_atom_index = ref_idx[0]
        else:
            self.reference_atom_index = -1
        self.constrained_indices = [constr.a for constr in self.mol_ase.constraints]
        self.natoms = len(self.mol_ase)


    def generate_molecule_rotations(
            self, 
            x_rot_angles : list,
            y_rot_angles : list,
            z_rot_angles : list,
            vert_angles_list : list = [0],
            individual_rotations : list = None,
            save_image = False,
            VERBOSE : bool = False
            ):
        '''
        Returns list (of ase Atoms) of all the rotated configurations and their labels.

        Args:
        - x_rot_angles: list of angles to rotate around the x axis
        - y_rot_angles: list of angles to rotate around the y axis
        - z_rot_angles: list of angles to rotate around the z axis
        - vert_angles_list: custom list of rotation angles for when the molecule is vertical (for y_angle = +/- 90°)
        - save_image: decide wether to save an image of all the molecular rotations
        '''
        
        if VERBOSE: print('Generating molecular configurations...')

        configs_ase = []
        labels = []


        #just to enter the loop and generate the (only one) unmodified configuration
        if not x_rot_angles : x_rot_angles = [0] #x_rot
        if not y_rot_angles : y_rot_angles = [0] #y_rot
        if not z_rot_angles : z_rot_angles = [0] #z_rot
        
        #build list of rotations (all combinations of x,y,z)
        rotations_list = []
        for y_angle in y_rot_angles:
            if (y_angle == 90 or y_angle == -90):
                for ang in vert_angles_list:
                    rotations_list.append([ang, y_angle, 0.0]) #x rotation is equivalent to z rotation              
                continue
            for x_angle in x_rot_angles:         
                for z_angle in z_rot_angles:
                        rotations_list.append([x_angle, y_angle, z_angle])

        #add also the individual rotations
        if individual_rotations:
            rotations_list += individual_rotations
                                
        #perform the rotations 
        for (x_angle, y_angle, z_angle) in rotations_list:
                
            mol = self.mol_ase.copy()
           
            mol.rotate(x_angle, 'x')
            mol.rotate(y_angle, '-y')
            mol.rotate(z_angle, 'z')          

            configs_ase.append(mol)
            labels.append( '{0},{1},{2},'.format(x_angle, y_angle, z_angle) )
                
        if VERBOSE: print('All molecular configurations generated.') 

        if(save_image):
            save_rotations_images(configs_ase, labels, VERBOSE=True)

        return configs_ase, labels


def save_rotations_images(configs_ase : list, labels : list, figname : str = 'molecule_orientations.png', VERBOSE : bool = False):
    '''
    Plot all the rotated molecules from top view

    Args:
    - configs_ase: list of ase Atoms objects containing all the rotated molecules
    - labels: list of strings containing the info on each rotation
    '''

    if VERBOSE: print("Saving image to {0}".format(figname))

    rows_fig = max(int(np.ceil(len(configs_ase)/5)), 1)
    cols_fig = max(int(np.ceil(len(configs_ase)/rows_fig)), 1)
    fig = plt.figure(figsize=(3*10 * cols_fig / 5, 3*5 * rows_fig / 3))
    axes = [fig.add_subplot(rows_fig,cols_fig,i) for i in range(1,len(configs_ase) + 1)]

    for i, atoms in enumerate(configs_ase):
        conf = atoms.copy()
        center = (conf.cell[:][0]/2 + conf.cell[:][1]/2)
        conf.translate(center)
        plot_atoms(conf, axes[i], show_unit_cell=2)
        axes[i].set_title('({0}, {1}, {2})'.format(* labels[i].split(',')))
    fig.suptitle('Molecule orientations (xrot, yrot, zrot)')
    fig.savefig(figname, dpi=800, bbox_inches='tight')

    if VERBOSE: print("Image saved.")