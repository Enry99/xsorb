#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue 28 Feb 2023
@author: Enrico Pedretti

Helper class to manage the molecule

"""

import numpy as np
from pymatgen.io import ase
from ase.io import read
from matplotlib import pyplot as plt
from ase.visualize.plot import plot_atoms


class Molecule:

    def __init__(self, molecule_filename : str, molecule_axis_atoms : list = None, axis_vector: list = None, atoms_subset : list = None):
        '''
        Reads molecule from file (e.g. Quantum ESPRESSO pwi/pwo or .xyz)
        Optional parameters:
        -atoms_subset: indices of the atoms to include. Only removes atom NOT in this list, does not check
        if the atoms in the list actually exist in the molecule
        '''

        self.mol_ase = read(molecule_filename)

        #align axis to x axis
        if molecule_axis_atoms:
            x1 = np.array(self.mol_ase.get_positions()[molecule_axis_atoms[0]])
            x2 = np.array(self.mol_ase.get_positions()[molecule_axis_atoms[1]])
            self.mol_ase.rotate(x2-x1, 'x')
        elif axis_vector:
            self.mol_ase.rotate(axis_vector, 'x')

        self.mol_pymat = ase.AseAtomsAdaptor.get_molecule(self.mol_ase, charge_spin_check=False)


        #select atoms subset
        self.original_mol_ase = self.mol_ase.copy() #before removing atoms
        all_indices = range(0, len(self.mol_ase))
        if atoms_subset:            
            remove_indices = [index for index in all_indices if index not in atoms_subset]
            self.reindex_map = [index for index in all_indices if index in atoms_subset]
            self.mol_pymat.remove_sites(remove_indices)
            self.mol_ase = ase.AseAtomsAdaptor.get_atoms(self.mol_pymat)
        else:
            self.reindex_map = [*all_indices]

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
            vert_rotations : list= [0],
            screw_rotations : list=[0],
            horiz_rotations : list=[0],
            no_vert_rotx : bool = False,
            distance_from_surf = 2., 
            min_distance=1., 
            save_image = False
            ):
        '''
        Returns a pymatgen list of all the rotated configurations and their labels.

        Args:
            atom_index, coords: index of the atom to be translated to the origin, which will be placed on the selected site. Alternatively you can specify the [x,y,z] coordinates that you want to place in the origin.
            horiz_rotations, vert_rotations: list of the rotations. vertical means around y axis (performed first), horizontal around z.
            min_distance: minimum distance that any atom can have from the surface. The whole molecule is translated upwards if necessary, in order to satisfy this condition.
        '''
        
        configs_ase = []
        labels = []   #format [[xrot, yrot, zrot], z_coord]
        delta_z_values = []

        transl_molecule = self._set_atom_to_origin(atom_index, coords)

        #just to enter the loop and generate the (only one) unmodified configuration
        if not vert_rotations : vert_rotations = [0]
        if not screw_rotations : screw_rotations = [0]
        if not horiz_rotations : horiz_rotations = [0]
        
        for vert_angle in vert_rotations:
                     
            j = -1 #index for the pre-relax distances (no horiz. rotations)
            for i, screw_angle in enumerate(screw_rotations):             

                j += 1
                if (isinstance(distance_from_surf, float)):
                    distance = distance_from_surf
                else: distance = distance_from_surf[j]                

                for hor_angle in horiz_rotations:

                    mol = transl_molecule.copy()

                    mol.rotate(screw_angle, 'x')
                    mol.rotate(vert_angle, '-y')
                    mol.rotate(hor_angle, 'z')

                    mol.translate([0,0,distance])

                    #Check for min_distance and translate accordingly
                    assert(min_distance <= distance)
                    atoms_too_close = [atom.z for atom in mol if atom.z < min_distance]
                    if atoms_too_close:
                        z_min = min(atoms_too_close)                       
                        delta_z = -z_min + min_distance
                        if(False): delta_z_values.append(delta_z)

                        #print('Config. roty={0} rotz={1}: {2} atoms closer than intended distance. Translating molecule {3:.3f} upwards '.format(vert_angle, hor_angle, len(atoms_too_close), delta_z))
                        mol.translate([0,0,delta_z])              

                    configs_ase.append(mol)

                    if(False): labels.append( ['{0},{1},{2},'.format(screw_angle, vert_angle, hor_angle), 0] )
                    labels.append( ['{0},{1},{2},'.format(screw_angle, vert_angle, hor_angle), '{:.3f}'.format(mol.get_positions()[atom_index][2])] )

                    if vert_angle == 90: break
                
                if vert_angle == 90 and no_vert_rotx: break

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

        #Conversion into pymatgen object
        configs_pymat = [ase.AseAtomsAdaptor.get_molecule(mol_ase, charge_spin_check=False) for mol_ase in configs_ase]

        return configs_pymat, labels
