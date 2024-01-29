#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: Enrico Pedretti

Small helper class to collect the settings read from settings.in

"""


import sys
import numpy as np
import input

#OK (code agnostic) 
class Settings:
    '''
    Class to read the settings.in and store all the parameters for the adsorption structures generation 
    and for the dft code

    Args:
    - settings_filename: path for the settings file
    - read_energies: try to read molecule and slab energies from files
    '''

    def __init__(self, settings_filename : str = "settings.in", read_energies : bool = True, VERBOSE : bool = True):

        #variables read from settings.in
        script_settings_dict, self.dftprogram_settings_dict = input.read_input_file(settings_filename)
    

        #check for existence of the main blocks
        blocks = ['INPUT','STRUCTURE']
        for block in blocks:
            if block not in script_settings_dict:
                raise RuntimeError(f"{block} block not found in settings.in.") 


        #check existence of mandatory flags in blocks
        mandatory_flags = {
            'INPUT' : ['slab_filename','molecule_filename', 'jobscript'], 
            'STRUCTURE' : ['selected_atom_index', 'molecule_axis', 'x_rot_angles', 'y_rot_angles', 'z_rot_angles'], 
            }
        for block in mandatory_flags:
            for flag in mandatory_flags[block]:
                if flag not in script_settings_dict[block]:
                    raise RuntimeError(f"Mandatory flag {flag} not found in block {block} while reading settings.in.")


        script_settings_dict['INPUT']['jobscript'] = [x.strip("'") for x in script_settings_dict['INPUT']['jobscript'].split()]
        script_settings_dict['STRUCTURE']['molecule_axis'] = [x.strip("'") for x in script_settings_dict['STRUCTURE']['molecule_axis'].split()]

        if script_settings_dict['STRUCTURE']['molecule_axis'][0] == 'atoms':
            if len(script_settings_dict['STRUCTURE']['molecule_axis'][1:]) != 2:
                print("Error: you must specify two integers as indices of the two atoms for the molecule axis.")
                sys.exit(1)
        elif script_settings_dict['STRUCTURE']['molecule_axis'][0] == 'vector':
            if len(script_settings_dict['STRUCTURE']['molecule_axis'][1:]) != 3:
                raise ValueError("Error: you must specify three numbers for the vector representing the molecule axis.")         
        else:
            raise ValueError("Error: molecule_axis must be specified with 'atoms' or 'vector'")


        if 'vertical_angles' in script_settings_dict['STRUCTURE']:
            script_settings_dict['STRUCTURE']['vertical_angles'] = [x.strip("'") for x in script_settings_dict['STRUCTURE']['vertical_angles'].split()]
            if (script_settings_dict['STRUCTURE']['vertical_angles'][0] != 'x' 
                and script_settings_dict['STRUCTURE']['vertical_angles'][0] != 'z'
                and script_settings_dict['STRUCTURE']['vertical_angles'][0] != 'none'
                and script_settings_dict['STRUCTURE']['vertical_angles'][0] != 'list'):
                    raise ValueError("Error: vertical_angles requires either 'x', 'z', 'none' or 'list' [list of custom angles].")         

            if (script_settings_dict['STRUCTURE']['vertical_angles'][0] == 'x' 
                or script_settings_dict['STRUCTURE']['vertical_angles'][0] == 'z'
                or script_settings_dict['STRUCTURE']['vertical_angles'][0] == 'none') and len(script_settings_dict['STRUCTURE']['vertical_angles']) > 1 :
                    raise ValueError("Error: the specified mode for vertical_angles does not support additional arguments.")
            
            if(script_settings_dict['STRUCTURE']['vertical_angles'][0] == 'list' and len(script_settings_dict['STRUCTURE']['vertical_angles']) == 1):
                    raise ValueError("Error: the specified mode for vertical_angles requires at least one angle (e.g. 0)")
        else: script_settings_dict['STRUCTURE']['vertical_angles'] = ['x']

        if 'fixed_indices_slab' in script_settings_dict['STRUCTURE'] and 'fixed_layers_slab' in script_settings_dict['STRUCTURE']:
            raise ValueError("You can specify the fixed slab atoms either by 'fixed_indices_slab' or 'fixed_layers_slab', not both.")


        #set non-specified flags to default.
        optional_structue_flags_list = {
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
            'mol_before_slab'          : False,
            'sort_atoms_by_z'          : True,
            'translate_slab'           : True
        }
        for flag in optional_structue_flags_list:
            if flag not in script_settings_dict['STRUCTURE']:
                script_settings_dict['STRUCTURE'].update({flag : optional_structue_flags_list[flag]})


        #Finally read data from dictionary and initialize the actual variables

        #&INPUT:
        self.slab_filename          = script_settings_dict['INPUT']['slab_filename']
        self.molecule_filename      = script_settings_dict['INPUT']['molecule_filename']
        self.jobscript              = script_settings_dict['INPUT']['jobscript'][0]
        self.sbatch_command         = script_settings_dict['INPUT']['jobscript'][1]
        if 'screening_conv_thr' in script_settings_dict['INPUT']:
            self.screening_conv_thr = np.array(script_settings_dict['INPUT']['screening_conv_thr'].split(), dtype=float).tolist()
        else: self.screening_conv_thr = None
        self.E_slab_mol = np.array(
            script_settings_dict['INPUT']['E_slab_mol'].split(), dtype=float
            ).tolist() if 'E_slab_mol' in script_settings_dict['INPUT'] else []

        if len(self.E_slab_mol) != 2 and len(self.E_slab_mol) != 0:
            print("If you specify the tag E_slab_mol you must provide TWO values. Quitting.")
            sys.exit(1)

        if(not self.E_slab_mol and read_energies):
            if script_settings_dict['STRUCTURE']['mol_subset_atoms'] == ' ':
                print('Using a subset of the molecule atoms. Molecule energy not extracted.')
            else:
                try:
                    from ase.io import read
                    import ase_custom
                    #try to read E_salb_mol directly from files
                    slab_en = read(self.slab_filename).get_potential_energy()
                    mol_en = read(self.molecule_filename).get_potential_energy()
                    self.E_slab_mol = [slab_en, mol_en]
                    print("Slab and molecule energies read from files.")
                except:
                    print('It was not possible to obtain slab and molecule energy in any way. Total energies will be displayed instead.')
                    self.E_slab_mol = [0, 0]

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

        self.sort_atoms_by_z        = script_settings_dict['STRUCTURE']['sort_atoms_by_z']
        self.translate_slab         = script_settings_dict['STRUCTURE']['translate_slab']
        self.mol_before_slab        = script_settings_dict['STRUCTURE']['mol_before_slab']

        self.sites_find_args = {
            'symm_reduce':self.symm_reduce, 
            'near_reduce':self.near_reduce, 
            'no_obtuse_hollow':True}

        #import json
        #print(json.dumps(self.__dict__, indent=4))

