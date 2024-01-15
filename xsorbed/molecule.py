#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue 28 Feb 2023
@author: Enrico Pedretti

Small helper class to manage the molecule

"""

import numpy as np
from ase.io import read
from matplotlib import pyplot as plt
from ase.visualize.plot import plot_atoms
from ase.constraints import FixCartesian


import ase_custom 

class Molecule:

    def __init__(self, molecule_filename : str, molecule_axis_atoms : list = None, axis_vector: list = None, atoms_subset : list = None, fixed_indices_mol : list = None, fix_mol_xyz : list = None):
        '''
        Reads molecule from file (e.g. Quantum ESPRESSO pwi/pwo or .xyz)
        Optional parameters:
        -atoms_subset: indices of the atoms to include. Only removes atom NOT in this list, does not check
        if the atoms in the list actually exist in the molecule
        '''

        print('Loading molecule...')
        self.mol_ase = read(molecule_filename, results_required=False)

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
            c = [FixCartesian(atom_index, mask=[not x for x in fix_mol_xyz]) for atom_index in fixed_indices_mol]  #we need to negate: in qe 0 = fix, here 1(true)=fix
            self.mol_ase.set_constraint(c)
        ###############################################################
            
        #select atoms subset
        self.original_mol_ase = self.mol_ase.copy() #before removing atoms
        all_indices = range(0, len(self.mol_ase))
        if atoms_subset:            
            remove_indices = [index for index in all_indices if index not in atoms_subset]
            self.reindex_map = [index for index in all_indices if index in atoms_subset]
            del self.mol_ase[remove_indices]
        else:
            self.reindex_map = [*all_indices]

        self.natoms = len(self.mol_ase)
        
        print('Molecule loaded.')

    #translation and rotation functions defs
    def _set_atom_to_origin(self, atom_index : int = None, coords : list = None):
        '''
        Returns a copy of the molecule with the selected atom translated to the origin,
        or alternatively translated by the vector specified by coords
        '''
        mol = self.mol_ase.copy()

        if atom_index is not None:
            #Reindex if only some atoms were selected from input
            atom_coords = self.original_mol_ase.get_positions()[atom_index].tolist()
            #print(atom_coords)   
            #print(self.mol_ase.get_positions().tolist())     
            atom_reindex = self.mol_ase.get_positions().tolist().index(atom_coords)        
            mol.translate(-mol.get_positions()[atom_reindex])
        else:
            mol.translate(-np.array(coords)) #NOTE!: for now NOT implemented
        return mol

    def generate_molecule_rotations(
            self, 
            atom_index : int=None, 
            coords : list=None, 
            y_rot_angles : list= [0.0],
            x_rot_angles : list=[0.0],
            z_rot_angles : list=[0.0],
            vert_angles_list : list = [0.0],
            distance_from_surf : float = 2., 
            min_distance : float = 1.5, 
            save_image = False,
            ):
        '''
        Returns list (of ase Atoms) of all the rotated configurations and their labels.
        '''
        
        print('Generating molecular configurations...')
        configs_ase = []
        labels = []   #format [[xrot, yrot, zrot], z_coord]
        delta_z_values = []

        transl_molecule = self._set_atom_to_origin(atom_index, coords)

        #just to enter the loop and generate the (only one) unmodified configuration
        if not y_rot_angles : y_rot_angles = [0]   #y_rot
        if not x_rot_angles : x_rot_angles = [0] #x_rot
        if not z_rot_angles : z_rot_angles = [0] #z_rot
        

        rotations_list = []

        for y_angle in y_rot_angles:
            if (y_angle == 90 or y_angle == -90):
                for ang in vert_angles_list:
                    rotations_list.append([ang, y_angle, 0.0])               
                continue
            for x_angle in x_rot_angles:         
                for z_angle in z_rot_angles:
                        rotations_list.append([x_angle, y_angle, z_angle])
                                
            
        for (x_angle, y_angle, z_angle) in rotations_list:
                
            mol = transl_molecule.copy()

            mol.rotate(x_angle, 'x')
            mol.rotate(y_angle, '-y')
            mol.rotate(z_angle, 'z')

            mol.translate([0,0,distance_from_surf])

            #Check for min_distance and translate accordingly
            assert(min_distance <= distance_from_surf)
            atoms_too_close = [atom.z for atom in mol if atom.z < min_distance]
            if atoms_too_close:
                z_min = min(atoms_too_close)                       
                delta_z = -z_min + min_distance
                if(False): delta_z_values.append(delta_z)

                #print('Config. roty={0} rotz={1}: {2} atoms closer than intended distance. Translating molecule {3:.3f} upwards '.format(y_angle, z_angle, len(atoms_too_close), delta_z))
                mol.translate([0,0,delta_z])              

            configs_ase.append(mol)

            if(False): labels.append( ['{0},{1},{2},'.format(x_angle, y_angle, z_angle), 0] )
            labels.append( ['{0},{1},{2},'.format(x_angle, y_angle, z_angle), '{:.3f}'.format(mol.get_positions()[atom_index][2])] )
                


        if(False): #translation for all same z
            for i, mol in enumerate(configs_ase):
                mol.translate([0,0, max(delta_z_values)])
                labels[i][1] = '{:.3f}'.format(selected_atom_distance + mol.get_positions()[atom_index][2])
                
    
        if(save_image): #Plot all the configs from top view
            rows_fig = max(int(np.ceil(len(configs_ase)/5)), 1)
            cols_fig = max(int(np.ceil(len(configs_ase)/rows_fig)), 1)
            fig = plt.figure(figsize=(10 * cols_fig / 5, 5 * rows_fig / 3))
            axes = [fig.add_subplot(rows_fig,cols_fig,i) for i in range(1,len(configs_ase) + 1)]

            for i, conf in enumerate(configs_ase):
                plot_atoms(conf, axes[i])
                axes[i].set_title('({0}, {1}, {2})'.format(* labels[i][0].split(',')))
            fig.suptitle('Molecule orientations (xrot, yrot, zrot)')
            fig.savefig('molecule_orientations.png', dpi=800, bbox_inches='tight')
        
        print('All molecular configurations generated.')

        return configs_ase, labels
