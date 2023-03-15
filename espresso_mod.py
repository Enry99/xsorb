#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue 28 Feb 2023
@author: Enrico Pedretti

Subclass Espresso to modify write pwi in order to fix atoms xy only

"""

from ase.calculators.espresso import Espresso
from ase.calculators.calculator import FileIOCalculator
from ase.io import write
import numpy as np

class Espresso_mod(Espresso):

    def set_fixed_atoms(self, slab_atoms_indices: list, slab_remapping: list, mol_atom_indices: list, mol_remapping: list, natoms_slab : int, natoms_mol : int, fix_slab_xyz : list = [0,0,0], fix_mol_xyz : list = [0,0,1]):
        '''
        Sets the indices of the atoms that need to be fixed during dynamics.
        Args: 
            slab_atoms_indices: list of the slab atoms indices to be fixed
            mol_atom_index: index of the atom of the mol to be fixed (index in the slab+mol cell)
            FixZ_slab, FixZ_molecule: fix also the z coordinate
        '''

        self.slab_atoms_indices = slab_atoms_indices.copy()
        self.slab_remapping     = slab_remapping.copy()
        self.mol_atom_indices   = mol_atom_indices.copy()
        self.mol_remapping      = mol_remapping.copy()
        self.natoms_slab        = natoms_slab
        self.natoms_mol         = natoms_mol
        self.fix_slab_xyz       = fix_slab_xyz.copy()
        self.fix_mol_xyz        = fix_mol_xyz.copy()


    def set_starting_magnetizations(self, starting_magnetizations : list):

        if(starting_magnetizations):
            self.starting_mag = starting_magnetizations.copy()


    def set_vdwd2_coeffs(self, vdw2_c6 : list, vdw2_rvdw : list):
        if(vdw2_c6):
            self.vdw2_c6_list = vdw2_c6.copy()
        if(vdw2_rvdw):
            self.vdw2_rvdw_list = vdw2_rvdw.copy()


    def _fix_atoms_pwi(self):
        #Edit pwi after FileIOCalculator.write_input, adding the constraints on atom dynamics
        
        with open(self.label + '.pwi', 'r') as file:
            data = file.readlines()

        first_atom_index = data.index('ATOMIC_POSITIONS angstrom\n') + 1


        if -1 in self.slab_atoms_indices:
            for i in range(self.natoms_slab):
                data[first_atom_index+i] = data[first_atom_index+i].split('\n')[0] + '    {0} {1} {2}\n'.format( *self.fix_slab_xyz)
        else:
            for i in self.slab_atoms_indices:
                data[first_atom_index+self.slab_remapping.tolist().index(i)] = data[first_atom_index+self.slab_remapping.tolist().index(i)].split('\n')[0] + '    {0} {1} {2}\n'.format( *self.fix_slab_xyz)
        

        if -1 in self.mol_atom_indices:
            for i in range(self.natoms_slab, self.natoms_slab + self.natoms_mol):
                data[first_atom_index+i] = data[first_atom_index+i].split('\n')[0] + '    {0} {1} {2}\n'.format( *self.fix_mol_xyz) 
        else:
            for i in self.mol_atom_indices:
                data[first_atom_index + self.natoms_slab + self.mol_remapping.index(i)] = data[first_atom_index+ self.natoms_slab +self.mol_remapping.index(i)].split('\n')[0] + '    {0} {1} {2}\n'.format( *self.fix_mol_xyz)            

        with open(self.label + '.pwi', 'w') as file:
            file.writelines( data )


    def _starting_magnetization_pwi(self):
        #Edit pwi after FileIOCalculator.write_input, adding the initial magnetization, since otherwise it is overwritten.
        
        with open(self.label + '.pwi', 'r') as file:
            data = file.readlines()

        for i, mag in enumerate(self.starting_mag):

            for j, line in enumerate(data):
                if 'starting_magnetization({0})'.format(i+1) in line:
                    line_index = j         
            
            if mag is not None:
                data[line_index] = '   starting_magnetization({0}) = {1}\n'.format(i+1, mag)

        with open(self.label + '.pwi', 'w') as file:
            file.writelines( data )


    def _vdw2_coeffs_pwi(self):

        with open(self.label + '.pwi', 'r') as file:
            data = file.readlines()

        line_index=0
        for i, line in enumerate(data):
            if 'vdw_corr' in line:
                line_index = i+1
        
        if line_index == 0: 
            print('No vdw option found. vdw2 parameters NOT written in pwi.')
            return   #no vdw option found

        for i, coeff in enumerate(self.vdw2_c6_list):
            if coeff is not None:
                data.insert(line_index, '   london_c6({0})     = {1}\n'.format(i+1, coeff))

        for i, coeff in enumerate(self.vdw2_rvdw_list):      
            if coeff is not None:
                data.insert(line_index, '   london_rvdw({0})   = {1}\n'.format(i+1, coeff))


        with open(self.label + '.pwi', 'w') as file:
            file.writelines( data )



    #Redefinition of function from Espresso to add option to fix atoms in pwi.
    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        write(self.label + '.pwi', atoms, **self.parameters)
        if(hasattr(self, 'slab_atoms_indices')): self._fix_atoms_pwi()
        if(hasattr(self, 'starting_mag')): self._starting_magnetization_pwi()
        if(hasattr(self, 'vdw2_c6_list') or hasattr(self, 'vdw2_rvdw_list')): self._vdw2_coeffs_pwi()