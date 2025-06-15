#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Created on Tue 28 Feb 2023
# @author: Enrico Pedretti

"""
Module that contains the class Molecule, used to generate the rotations
"""

from __future__ import annotations
import warnings

import numpy as np
from ase.constraints import FixCartesian
from ase.build import connected_indices, split_bond

from xsorb.ase_custom.atoms import AtomsCustom
from xsorb.visualize.plot import plot_rotations_images
from xsorb.adsorptiondata.adsorptionstructure import MoleculeRotation, SurroundingSite


class Molecule:
    '''
        Class to generate the rotations of the molecule around a reference atom.
        The molecule is aligned to the x-axis, and the rotations are performed around
        the reference atom, which is translated to the origin.
        The molecule can be rotated around the x, y and z axes, with the angles specified
        in the `generate_molecule_rotations` method.

        Initialization parameters:
        - mol: Atoms object of the molecule
        - atom_indexes: list of indices of the reference atom (indexing starting from 0,
            in the order of the input file). If -1, use geometry center of the molecule.
        - molecule_axis_mode: 'atoms' or 'vector'
        - molecule_axis_values: values defining the x-axis of the molecule:
            'atoms' mode: indices of the two atoms (a = r2 - r1)
            'vector' mode: [ax, ay, az] vector that will be aligned to the x-axis.
        - atoms_subset: indices of the atoms to include. Only removes atom NOT in this list,
            does not check if the atoms in the list actually exist in the molecule
        - break_bond_indices: indices of the two atoms of the bond to be broken.
            The fragment containing the first atom will be kept.
        - fixed_indices_mol: list of specific atoms to be fixed (-1 : fix all)
            (indices start from 0, with the ordering of the atoms in the input file)
        - fix_mol_xyz: which coordinates to fix for the fixed atoms,
            e.g. [True, True, False] = fix motion in x,y, free to move along z.
    '''


    def __init__(self, mol: AtomsCustom,
                 atom_indexes: list[int],
                 molecule_axis_mode: str,
                 molecule_axis_values : list,
                 atoms_subset : list[int] | None = None,
                 break_bond_indices : list[int] | None = None,
                 fixed_indices_mol : list[int] | None = None,
                 fix_mol_xyz : list[bool] | None = None):


        self.mol_ase = mol.copy()

        #align axis to x axis
        if molecule_axis_mode == 'atoms':
            if len(molecule_axis_values) != 2:
                raise ValueError('molecule_axis_values must contain exactly two atom indices.')
            x1 = self.mol_ase.positions[molecule_axis_values[0]]
            x2 = self.mol_ase.positions[molecule_axis_values[1]]
            self.mol_ase.rotate(x2-x1, 'x')
        elif molecule_axis_mode == 'vector':
            if len(molecule_axis_values) != 3:
                raise ValueError('molecule_axis_values must contain exactly three values.')
            self.mol_ase.rotate(molecule_axis_values, 'x')
        else:
            raise ValueError('molecule_axis_mode must be either "atoms" or "vector".')

        #SET CONSTRAINTS
        if fixed_indices_mol:
            if -1 in fixed_indices_mol:
                fixed_indices_mol = list(range(len(self.mol_ase)))
            c = [FixCartesian(idx, mask=fix_mol_xyz) for idx in fixed_indices_mol]
                #True = fixed, False = free
            self.mol_ase.set_constraint(c)
        else: self.mol_ase.set_constraint() #clean possible constraints read from file
        ###############################################################

        #select atoms subset
        atom_indexes_positions = [self.mol_ase.positions[idx] \
                                  if idx != -1 else None for idx in atom_indexes]
        if atoms_subset:
            self.mol_ase = self.mol_ase[atoms_subset]

        elif break_bond_indices:
            if len(break_bond_indices) != 2:
                raise ValueError('break_bond_indices must contain exactly two atom indices.')

            # check if the two indices are connected
            connected = connected_indices(self.mol_ase, break_bond_indices[0])
            if break_bond_indices[1] not in connected:
                raise ValueError('The two atoms in break_bond_indices are not bonded.')

            # split the bond and keep the fragment containing the first atom
            self.mol_ase = split_bond(self.mol_ase, break_bond_indices[0],
                                      break_bond_indices[1])[0]

        #retrieve the correct index of the reference atoms after subsetting
        for i, atom_pos in enumerate(atom_indexes_positions):
            if atom_pos is not None:
                atom_indexes[i] = [atom.index for atom in self.mol_ase \
                                   if np.allclose(atom.position, atom_pos)][0]
            else:
                atom_indexes[i] = -1

        self.reference_atom_indices = atom_indexes

        #Translate reference atom to origin
        self.translated_mols = {}
        for atom_index in self.reference_atom_indices:
            mol = self.mol_ase.copy()
            if atom_index != -1:
                mol.translate(-self.mol_ase.positions[atom_index])
            else:
                #center of positions (not of mass)
                mol.translate(-self.mol_ase.positions.mean(axis=0))
            self.translated_mols.update({atom_index : mol})

        #set other useful attributes
        self.constrained_indices = [constr.index[0] for constr in self.mol_ase.constraints]
        self.natoms = len(self.mol_ase)


    def generate_molecule_rotations(self,
            which_index : int,
            x_rot_angles : list[float],
            y_rot_angles : list[float],
            z_rot_angles : list[float|SurroundingSite],
            vert_angles_list : list[float] | None,
            individual_rotations : list[float] | None = None,
            surrounding_exclude_main : bool = False,
            save_image : bool = False,
            verbose : bool = False
            ):
        '''
        Returns list (of ase Atoms) of all the rotated configurations and their labels.

        Args:
        - which_index: index of the reference atom to use as pivot for the rotations
        - x_rot_angles: list of angles to rotate around the x axis
        - y_rot_angles: list of angles to rotate around the y axis
        - z_rot_angles: list of angles to rotate around the z axis
        - vert_angles_list: custom list of rotation angles for when the molecule is vertical
            (for y_angle = +/- 90°)
        - save_image: decide wether to save an image of all the molecular rotations
        '''

        if verbose:
            print('Generating molecular configurations...')

        if which_index not in self.reference_atom_indices:
            raise ValueError('The index of the reference atom is not valid.')

        mol_rotations_ase : list[MoleculeRotation] = []


        #just to enter the loop and generate the (only one) unmodified configuration,
        # useful if the argument was []
        if not x_rot_angles:
            x_rot_angles = [0.] #x_rot
        if not y_rot_angles:
            y_rot_angles = [0.] #y_rot
        if not z_rot_angles:
            z_rot_angles = [0.] #z_rot

        #build list of rotations (all combinations of x,y,z)
        rotations_list : list[list[float|SurroundingSite]] = []
        for y_angle in y_rot_angles:
            if (int(y_angle) in [90, -90]): #special case, no further loops with x and z angles
                # z rotation is equivalent to x rotation
                rotations_list.extend([[ang, y_angle, 0] for ang in vert_angles_list] \
                                      if vert_angles_list else [[0, y_angle, 0]])
                continue
            for x_angle in x_rot_angles:
                for z_angle in z_rot_angles:
                    rotations_list.append([x_angle, y_angle, z_angle])

        #add also the individual rotations
        if individual_rotations:
            rotations_list += individual_rotations

        #perform the rotations
        for (x_angle, y_angle, z_angle) in rotations_list:

            mol = self.translated_mols[which_index].copy()

            mol.rotate(x_angle, 'x')
            mol.rotate(y_angle, '-y')

            if isinstance(z_angle, SurroundingSite):
                if z_angle.duplicate_main and surrounding_exclude_main:
                    continue
                vector = z_angle.vector.copy()
                vector[2] = 0 #project on xy plane
                mol.rotate('x', vector)
            else:
                mol.rotate(z_angle, 'z')

            mol_rotations_ase.append(MoleculeRotation(mol,
                                                      str(x_angle),
                                                      str(y_angle),
                                                      str(z_angle),
                                                      mol_atom=which_index))

        if verbose:
            print('All molecular configurations generated.')

        if save_image:
            plot_rotations_images(mol_rotations_ase, verbose=True)

        return mol_rotations_ase
