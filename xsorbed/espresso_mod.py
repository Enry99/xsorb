#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue 28 Feb 2023
@author: Enrico Pedretti

Subclass of class Espresso to modify pwi written by ASE

"""

from ase.calculators.espresso import Espresso
from ase.calculators.calculator import FileIOCalculator
from ase.io import write
import sys
import numpy as np

class Espresso_mod(Espresso):

    def __init__(self, filename='input.pwi', **kwargs):
        super().__init__(**kwargs)
        self.filename = filename
        self.pseudo = kwargs['pseudopotentials'].keys()


    def set_fixed_atoms(self, slab_atoms_indices: list, slab_remapping: list, mol_atom_indices: list, mol_remapping: list, natoms_slab : int, natoms_mol : int, fix_slab_xyz : list = [0,0,0], fix_mol_xyz : list = [0,0,1]):
        '''
        Sets the indices of the atoms that need to be fixed during dynamics.
        '''

        self.slab_atoms_indices = slab_atoms_indices.copy()
        self.slab_remapping     = slab_remapping.copy()
        self.mol_atom_indices   = mol_atom_indices.copy()
        self.mol_remapping      = mol_remapping.copy()
        self.natoms_slab        = natoms_slab
        self.natoms_mol         = natoms_mol
        self.fix_slab_xyz       = fix_slab_xyz.copy()
        self.fix_mol_xyz        = fix_mol_xyz.copy()


    def set_system_flags(self, starting_magnetizations : list, flags_i : list):
        self.starting_mag = starting_magnetizations.copy()
        self.flags_i = flags_i.copy()


    def _fix_atoms_pwi(self):
        #Edit pwi after FileIOCalculator.write_input, adding the constraints on atom dynamics
        
        with open(self.filename, 'r') as file:
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

        with open(self.filename, 'w') as file:
            file.writelines( data )


    def _system_flags_pwi(self):
        #Edit pwi after FileIOCalculator.write_input, adding the initial magnetization, since otherwise it is overwritten.
        
        with open(self.filename, 'r') as file:
            data = file.readlines()


        ASE_DID_WRITE_MAGNETIZATIONS = False
        for i, mag in enumerate(self.starting_mag):

            found = False
            for j, line in enumerate(data):
                if 'starting_magnetization({0})'.format(i+1) in line:
                    line_index = j
                    found = True
                    ASE_DID_WRITE_MAGNETIZATIONS = True     
            
            if found and mag is not None:
                data[line_index] = '   starting_magnetization({0}) = {1}\n'.format(i+1, mag)

        in_SYSTEM = False
        j0 = -1
        for j, line in enumerate(data):
            if 'SYSTEM' in line: in_SYSTEM = True
            if('/' in line and in_SYSTEM):
                j0 = j
                break

        if not ASE_DID_WRITE_MAGNETIZATIONS:
            for i, mag in enumerate(self.starting_mag):
                if mag is not None:
                    data.insert(j0, '   starting_magnetization({0}) = {1}\n'.format(i+1, mag)) 
                    j0 += 1

        for i, flag in enumerate(self.flags_i):
            data.insert(j0, '   {0} = {1}\n'.format(flag[0], flag[1]))
            j0 += 1


        
        with open(self.filename, 'w') as file:
            file.writelines( data )

    def _sort_pseudo(self):
        with open(self.filename, 'r') as file:
            data = file.readlines()

            for i, line in enumerate(data):
                if 'ATOMIC_SPECIES' in line.upper():
                    block = data[i+1 : i+len(self.pseudo)+1]
                    i0 = i
                    break

            if(np.any(np.array(block) == '\n')):
                print('Number of ATOMIC_SPECIES not compatible with input structures. Quitting.')
                sys.exit(1)
            
            for j, ps in enumerate(self.pseudo):
                data[i0+j+1] = block[[x.split()[0] for x in block].index(ps)]

        with open(self.filename, 'w') as file:
            file.writelines( data )


    #Redefinition of function from Espresso to add option to fix atoms in pwi.
    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        write(self.filename, atoms, **self.parameters)
        self._sort_pseudo()
        if(hasattr(self, 'slab_atoms_indices')): self._fix_atoms_pwi()
        if(hasattr(self, 'starting_mag')): self._system_flags_pwi()
