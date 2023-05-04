#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: Enrico Pedretti

Small helper class to collect the settings read from settings.in

"""


import sys
import numpy as np
import input

class Settings:

    def __init__(self, settings_filename : str = "settings.in", read_energies = True) -> None:

        #variables read from settings.in
        script_settings_dict, self.espresso_settings_dict, self.atomic_species, self.kpoints, self.last_dump = input.read_input_file(settings_filename)
    

        #check for existence of the main blocks
        blocks = ['INPUT','STRUCTURE']
        for block in blocks:
            if block not in script_settings_dict:
                print("'"+block+"' block not found in settings.in.")
                sys.exit(1)  


        #check existence of mandatory flags in blocks
        mandatory_flags = {
            'INPUT' : ['slab_filename','molecule_filename', 'jobscript'], 
            'STRUCTURE' : ['selected_atom_index', 'molecule_axis', 'x_rot_angles', 'y_rot_angles', 'z_rot_angles'], 
            }
        for block in mandatory_flags:
            for flag in mandatory_flags[block]:
                if flag not in script_settings_dict[block]:
                    print("Mandatory flag {0} not found in block {1} while reading settings.in.".format(flag, block))
                    sys.exit(1)

        script_settings_dict['INPUT']['jobscript'] = [x.strip("'") for x in script_settings_dict['INPUT']['jobscript'].split()]
        script_settings_dict['STRUCTURE']['molecule_axis'] = [x.strip("'") for x in script_settings_dict['STRUCTURE']['molecule_axis'].split()]

        if script_settings_dict['STRUCTURE']['molecule_axis'][0] == 'atoms':
            if len(script_settings_dict['STRUCTURE']['molecule_axis'][1:]) != 2:
                print("Error: you must specify two integers as indices of the two atoms for the molecule axis.")
                sys.exit(1)
        elif script_settings_dict['STRUCTURE']['molecule_axis'][0] == 'vector':
            if len(script_settings_dict['STRUCTURE']['molecule_axis'][1:]) != 3:
                print("Error: you must specify three numbers for the vector representing the molecule axis.")
                sys.exit(1)            
        else:
            print("Error: molecule_axis must be specified with 'atoms' or 'vector'")
            sys.exit(1)

        print(script_settings_dict['STRUCTURE']['vertical_angles'])
        if 'vertical_angles' in script_settings_dict['STRUCTURE']:
            script_settings_dict['STRUCTURE']['vertical_angles'] = [x.strip("'") for x in script_settings_dict['STRUCTURE']['vertical_angles'].split()]
            if (script_settings_dict['STRUCTURE']['vertical_angles'][0] != 'x' 
                and script_settings_dict['STRUCTURE']['vertical_angles'][0] != 'z'
                and script_settings_dict['STRUCTURE']['vertical_angles'][0] != 'none'
                and script_settings_dict['STRUCTURE']['vertical_angles'][0] != 'list'):
                    print("Error: vertical_angles requires either 'x', 'z', 'none' or 'list' [list of custom angles].")
                    sys.exit(1)             

            if (script_settings_dict['STRUCTURE']['vertical_angles'][0] == 'x' 
                or script_settings_dict['STRUCTURE']['vertical_angles'][0] == 'z'
                or script_settings_dict['STRUCTURE']['vertical_angles'][0] == 'none') and len(script_settings_dict['STRUCTURE']['vertical_angles']) > 1 :
                    print("Error: the specified mode for vertical_angles does not support additional arguments.")
                    sys.exit(1) 
            if(script_settings_dict['STRUCTURE']['vertical_angles'][0] == 'list' and len(script_settings_dict['STRUCTURE']['vertical_angles']) == 1):
                    print("Error: the specified mode for vertical_angles requires at least one angle (e.g. 0)")
                    sys.exit(1)
        else: script_settings_dict['STRUCTURE']['vertical_angles'] = ['x']

        if 'fixed_indices_slab' in script_settings_dict['STRUCTURE'] and 'fixed_layers_slab' in script_settings_dict['STRUCTURE']:
            print("You can specify the fixed slab atoms either by 'fixed_indices_slab' or 'fixed_layers_slab', not both.")
            sys.exit(1)


        #set non-specified flags to default.
        optional_flags_list = {
            'symm_reduce'              : 0.01,
            'near_reduce'              : 0.01,
            'surface_height'           : 0.9,
            'layers_height'            : 0.5,
            'selected_sites'           : ' ',
            'mol_subset_atoms'         : ' ',
            'screening_atom_distance'  : 2.0,
            'screening_min_distance'   : 1.5,
            'relax_atom_distance'      : 2.0,
            'relax_min_distance'       : 1.5,
            'fixed_indices_slab'       : ' ',
            'fixed_layers_slab'        : ' ',
            'fixed_indices_mol'        : ' ',
            'fix_slab_xyz'             : '0 0 0',
            'fix_mol_xyz'              : '0 0 1',
        }
        for flag in optional_flags_list:
            if flag not in script_settings_dict['STRUCTURE']:
                script_settings_dict['STRUCTURE'].update({flag : optional_flags_list[flag]})
        if 'ion_dynamics' not in self.espresso_settings_dict['IONS']:
            self.espresso_settings_dict['IONS'].update({'ion_dynamics': 'bfgs'})
        if len(self.kpoints) == 2: #so no gamma but also no offset
            self.kpoints.append([0,0,0])


        #Finally read data from dictionary and initialize the actual variables

        #&INPUT:
        self.slab_filename          = script_settings_dict['INPUT']['slab_filename']
        self.molecule_filename      = script_settings_dict['INPUT']['molecule_filename']
        self.jobscript              = script_settings_dict['INPUT']['jobscript'][0]
        self.sbatch_command         = script_settings_dict['INPUT']['jobscript'][1]
        self.E_slab_mol = np.array(
            script_settings_dict['INPUT']['E_slab_mol'].split(), dtype=float
            ).tolist() if 'E_slab_mol' in script_settings_dict['INPUT'] else []

        if len(self.E_slab_mol) != 2 and len(self.E_slab_mol) != 0:
            print("If you specify the tag E_slab_mol you must provide TWO values. Quitting.")
            sys.exit(1)

        if(not self.E_slab_mol and read_energies):
            #try to read E_salb_mol directly from files (if files are completed pwos)
            if(self.slab_filename.split('.')[-1] == 'pwo' 
                and self.molecule_filename.split('.')[-1] == 'pwo'
                and script_settings_dict['STRUCTURE']['mol_subset_atoms'] == ' '):

                slab_en = self._read_energy(self.slab_filename)
                mol_en  = self._read_energy(self.molecule_filename)
                if slab_en is not None and mol_en is not None:
                    print("Slab and molecule energies read from files.")
                    self.E_slab_mol = [slab_en, mol_en]


        #&STRUCTURE
        self.symm_reduce            = float(script_settings_dict['STRUCTURE']['symm_reduce'])
        self.near_reduce            = float(script_settings_dict['STRUCTURE']['near_reduce'])
        self.surface_height         = float(script_settings_dict['STRUCTURE']['surface_height'])
        self.layers_height          = float(script_settings_dict['STRUCTURE']['layers_height'])
        self.selected_sites         = np.array(script_settings_dict['STRUCTURE']['selected_sites'].split(), dtype=int).tolist()

        self.mol_subset_atoms       = np.array(script_settings_dict['STRUCTURE']['mol_subset_atoms'].split(), dtype=int).tolist()
        self.molecule_axis_atoms    = np.array(script_settings_dict['STRUCTURE']['molecule_axis'][1:], dtype=int
                                               ).tolist() if script_settings_dict['STRUCTURE']['molecule_axis'][0] == 'atoms' else []
        self.axis_vector            = np.array(script_settings_dict['STRUCTURE']['molecule_axis'][1:], dtype=float
                                               ).tolist() if script_settings_dict['STRUCTURE']['molecule_axis'][0] == 'vector' else []
        self.selected_atom_index    = int(script_settings_dict['STRUCTURE']['selected_atom_index'])
        self.screening_atom_distance= float(script_settings_dict['STRUCTURE']['screening_atom_distance'])
        self.screening_min_distance = float(script_settings_dict['STRUCTURE']['screening_min_distance'])
        self.relax_atom_distance    = float(script_settings_dict['STRUCTURE']['relax_atom_distance'])
        self.relax_min_distance     = float(script_settings_dict['STRUCTURE']['relax_min_distance'])                
        self.x_rot_angles           = np.array(script_settings_dict['STRUCTURE']['x_rot_angles'].split(), dtype=float).tolist()
        self.y_rot_angles           = np.array(script_settings_dict['STRUCTURE']['y_rot_angles'].split(), dtype=float).tolist()  
        self.z_rot_angles           = np.array(script_settings_dict['STRUCTURE']['z_rot_angles'].split(), dtype=float).tolist()

        if script_settings_dict['STRUCTURE']['vertical_angles'][0] == 'x':            
            self.vertical_angles    = self.x_rot_angles
        elif script_settings_dict['STRUCTURE']['vertical_angles'][0] == 'z':
            self.vertical_angles    = self.z_rot_angles
        elif script_settings_dict['STRUCTURE']['vertical_angles'][0] == 'none':
            self.vertical_angles    = np.array([0.])
        else:
            self.vertical_angles    = np.array(script_settings_dict['STRUCTURE']['vertical_angles'][1:], dtype=float).tolist()
         

        self.fixed_indices_slab     = np.array(script_settings_dict['STRUCTURE']['fixed_indices_slab'].split(), dtype=int).tolist()
        self.fixed_layers_slab      = np.array(script_settings_dict['STRUCTURE']['fixed_layers_slab'].split(), dtype=int).tolist()
        self.fixed_indices_mol      = np.array(script_settings_dict['STRUCTURE']['fixed_indices_mol'].split(), dtype=int).tolist()
        self.fix_slab_xyz           = np.array(script_settings_dict['STRUCTURE']['fix_slab_xyz'].split(), dtype=int).tolist() 
        self.fix_mol_xyz            = np.array(script_settings_dict['STRUCTURE']['fix_mol_xyz'].split(), dtype=int).tolist()


        #&PSEUDO
        self.pseudopotentials       = {}
        for i, typ in enumerate(self.atomic_species):
            self.pseudopotentials.update({typ[0] : typ[2]})

        #KPOINTS
        if(len(self.kpoints) > 1): #so no gamma
            self.kpoints[1] = np.array(self.kpoints[1], dtype=int).tolist()
            self.kpoints[2] = np.array(self.kpoints[2], dtype=int).tolist()


        #Espresso &CONTROL
        self.espresso_settings_dict['CONTROL'].update({'outdir' : 'WORK'}) #TODO: lasciarlo scegliere all'utente

        #Espresso &SYSTEM (correction to bugs in ase function to write pwis)
        self.starting_mag = [None]*len(self.pseudopotentials) #this needs special treatment because it must be overwritten, while the others are just appended
        self.flags_i = []
        keys_to_be_removed = []

        for key in self.espresso_settings_dict['SYSTEM']:          
            if '(' in key and ')' in key:
                i = int(key.split('(')[1].split(')')[0])
                if i > len(self.pseudopotentials):
                    print("Flag {0} out of range (> than the nubmer of specified ATOMIC_SPECIES). Quitting.".format(key))
                    sys.exit(1)

                if 'starting_magnetization' in key:               
                    self.starting_mag[i-1] = self.espresso_settings_dict['SYSTEM'][key]
                else:
                    self.flags_i.append([key, self.espresso_settings_dict['SYSTEM'][key]])
                 
                keys_to_be_removed.append(key)              

        for key in keys_to_be_removed:
            self.espresso_settings_dict['SYSTEM'].pop(key)

        #Espresso &IONS
        self.ion_dynamics = self.espresso_settings_dict['IONS']['ion_dynamics']

        self.sites_find_args = {
            "distance":0,   #the distance isn't set from here, but inside the molecule, by translating the mol upwards by the chosen amount
            'symm_reduce':self.symm_reduce, 
            'near_reduce':self.near_reduce, 
            'no_obtuse_hollow':True}

        #import json
        #print(json.dumps(self.__dict__, indent=4))
        print(self.vertical_angles)

    def _read_energy(self, filename : str):
        with open(filename, 'r') as f:
            pwo = f.readlines()

        toten = 0
        scf_terminated = False
        for line in pwo: #make sure to get the last one (useful in relaxations)
            if '!' in line: 
                toten = line.split()[4]
            if 'End of self-consistent calculation' in line:
                scf_terminated = True
        if(scf_terminated and toten != 0): 
            return float(toten)
        else: return None

