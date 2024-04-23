#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu 11 May 2023

@author: Enrico Pedretti

Generation of molecular fragments and dissociation energies calculation

"""

import json, os, shutil, copy
import numpy as np
import pandas as pd
from xsorbed.settings import Settings
from xsorbed.molecule import Molecule
from xsorbed.calculations import launch_screening, final_relax
from xsorbed.io_utils import get_calculations_results, restart_jobs
from xsorbed.common_definitions import *
from xsorbed.dftcode_specific import FRAGMENTS_OUT_FILE_PATHS, FRAGMENTS_IN_FILE_PATHS, SBATCH_POSTFIX_FRAGS, \
    override_settings_isolated_fragment, override_settings_adsorbed_fragment, Calculator
from xsorbed import ase_custom


TEST = False   #set to true for testing: prints sbatch command instead of actually launching jobs for the isolated fragments


def generate_isolated_fragment(settings : Settings, fragment_dict : dict):
    """
    Generate isolated fragment from a molecule based on the provided settings and a dictionary for the fragment.

    Args:
       - settings (Settings): The settings object containing the required parameters for fragment generation.
       - fragment_dict (dict): A dictionary containing information about the fragment to be generated.

    Returns:
        Molecule object representing the isolated fragment.
    """

        
    if "mol_subset_atoms" in fragment_dict and "break_bond_indices" in fragment_dict:
        raise ValueError("You can specify a fragment either by its atoms or by the broken bonds, not both.")
    if "mol_subset_atoms" not in fragment_dict and "break_bond_indices" not in fragment_dict:
        raise RuntimeError("You must specify the fragment either through mol_subset_atoms or break_bond_indices.")
    
            
    axis = [x.strip("'") for x in fragment_dict['molecule_axis'].split()]
    if axis[0] == 'atoms':
        if len(axis[1:]) != 2:
            raise ValueError("Error: you must specify two integers as indices of the two atoms for the molecule axis.")
    elif axis[0] == 'vector':
        if len(axis[1:]) != 3:
            raise ValueError("Error: you must specify three numbers for the vector representing the molecule axis.")

    mol = Molecule(molecule_filename=settings.molecule_filename,
                            atom_index=fragment_dict["selected_atom_index"],
                            molecule_axis_atoms=np.array(axis[1:], dtype=int).tolist() if axis[0] == 'atoms' else None,
                            axis_vector=np.array(axis[1:], dtype=float).tolist() if axis[0] == 'vector' else None, 
                            atoms_subset=fragment_dict["mol_subset_atoms"] if "mol_subset_atoms" in fragment_dict else None,
                            break_bond_indices=fragment_dict["break_bond_indices"] if "break_bond_indices" in fragment_dict else None,
                            fixed_indices_mol=fragment_dict["fixed_indices_mol"] if "fixed_indices_mol" in fragment_dict else None,
                            fix_mol_xyz=fragment_dict["fix_mol_xyz"] if "fix_mol_xyz" in fragment_dict else None)
    
    if not mol.mol_ase.cell: #when reading from xyz, cell might not be set
        mol.mol_ase.set_cell([30,30,30])
        mol.mol_ase.set_pbc([True,True,True])

    return mol
    

def generate_all_isolated_fragments(settings : Settings, fragments_dict : dict, VERBOSE : bool = True):
    
    structures = []
    for fragment_name, fragment_dict in fragments_dict["fragments"].items():

        if fragment_name == "mol":
            raise ValueError('Error: fragment name "mol" not allowed (reserved for whole molecule).')
        
        structures.append(generate_isolated_fragment(settings, fragment_dict))

    return structures


def write_fragments_inputs(settings : Settings, 
                           fragments_dict : dict, 
                           molecules : list,
                           OVERRIDE_SETTINGS : bool = True, 
                           INTERACTIVE : bool = True):
    """
    Write input files for fragments based on the given settings, fragments dictionary, and structures.

    Args:
       - settings (Settings): The settings object containing the calculation parameters.
       - fragments_dict (dict): A dictionary containing information about the fragments.
       - molecules (list): A list of Molecule objects corresponding to the fragments.
       - OVERRIDE_SETTINGS (bool, optional): Whether to override the settings for each fragment. Defaults to True.
       - INTERACTIVE (bool, optional): Whether to interactively prompt the user for input. Defaults to True.

    Returns:
        list: A list of fragment names for which input files were written.
    """    

    os.makedirs('fragments', exist_ok=True)
    main_dir = os.getcwd()
    
    frag_list = []
    for fragment_name, molecule in zip(fragments_dict["fragments"], molecules):

        print(f"Writing input for fragment {fragment_name}...")
        outdir = f'fragments/{fragment_name}'
        os.makedirs(outdir, exist_ok=True)

        if INTERACTIVE:
            if os.path.exists(FRAGMENTS_OUT_FILE_PATHS[settings.program].format(fragment_name)): 
                print(f'{fragment_name} output file already present, possibly from a running calculation. You can decide to re-calculate it or skip it.')        
                while True:
                    answer = input('Re-calculate? ("y" = yes, "n" = no): ')
                    if answer == 'yes' or answer == 'y' or answer == 'no' or answer == 'n': 
                        break
                    else: print('Value not recognized. Try again.')
                if 'n' in answer: continue

            if os.path.exists(FRAGMENTS_IN_FILE_PATHS[settings.program].format(fragment_name)):
                print(f'{fragment_name} input file already present.')
                while True:
                    answer = input('Overwrite? ("y" = yes, "n" = no): ')
                    if answer == 'yes' or answer == 'y' or answer == 'no' or answer == 'n': 
                        break
                    else: print('Value not recognized. Try again.')
                    
                if 'n' in answer: continue

        if "dft_settings_override" in fragments_dict["fragments"][fragment_name]:
            if "vacuum" in fragments_dict["fragments"][fragment_name]["dft_settings_override"]:
                manual_dft_override = fragments_dict["fragments"][fragment_name]["dft_settings_override"]["vacuum"]
        else: manual_dft_override = None

        
        newsettings = copy.deepcopy(settings)
        if OVERRIDE_SETTINGS: 
            override_settings_isolated_fragment(newsettings, molecule.natoms, manual_dft_override)
        
        os.chdir(outdir)
        calc = Calculator(newsettings, fragment_name, molecule.mol_ase, '.') 
        calc.write_input(molecule.mol_ase)
        os.chdir(main_dir)

        frag_list.append(fragment_name)

    if INTERACTIVE: print('All input files written.') 

    return frag_list


def launch_fragments_jobs(program : str, jobscript : str, sbatch_command : str, fragments_list : list):
    '''
    Launch the calculations.

    Args:
    - program: DFT program. Possible values: 'ESPRESSO' or 'VASP'
    - jobscript: path of the jobscript file
    - sbatch_command: command to submit the jobscript (in Slurm it is sbatch)
    - fragments_list: list of the names of the fragments to launch
    '''
    main_dir = os.getcwd()

    for fragment in fragments_list:

        j_dir = f'fragments/{fragment}'
        os.makedirs(j_dir, exist_ok=True) #for QE (for VASP they are already created by write_inputs)
        shutil.copyfile(jobscript, f'{j_dir}/{jobscript_stdname}')
        
        os.chdir(j_dir)   ####################
        
        #change job title (only for slumr jobscripts)
        with open(jobscript_stdname, 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if "job-name" in line:
                    lines[i] = f"{line.split('=')[0]}={fragment}\n"
                    break
        with open(jobscript_stdname, 'w') as f:       
            f.writelines(lines)


        launch_string = f"{sbatch_command} {jobscript_stdname} {SBATCH_POSTFIX_FRAGS[program].format(fragment)}"
        if(TEST): print(launch_string)
        else: os.system(launch_string)  #launches the jobscript in j_dir from j_dir
        os.chdir(main_dir) ####################


def isolated_fragments(RUN=False):
    """
    Generate isolated fragments and launch fragment jobs.

    Args:
        - RUN (bool): Flag indicating whether to run the generated fragment jobs.

    """
    with open(framgents_filename, "r") as f:
        fragments_dict = json.load(f)

    settings = Settings()

    molecules = generate_all_isolated_fragments(settings, fragments_dict, VERBOSE=True)

    frag_list = write_fragments_inputs(settings, fragments_dict, molecules, OVERRIDE_SETTINGS=True, INTERACTIVE=True)

    if RUN: launch_fragments_jobs(settings.program, settings.jobscript, settings.sbatch_command, frag_list)




def edit_fragment_settings(lines : list, 
                           settings_dict : dict,
                           program : str):
    
    settings_lines_frag = [] #where to put all the lines

    script_section = ['@SETTINGS', '@/SETTINGS']
    dft_section = [f'@{program}', f'@/{program}']

    for i,line in enumerate(lines):
        if script_section[0] in line: 
            i_script_beg = i
        if script_section[1] in line:
            i_script_end = i
            break
    script_lines = lines[i_script_beg:i_script_end+1]

    for i,line in enumerate(lines):
        if dft_section[0] in line: 
            i_dft_beg = i
        if dft_section[1] in line:
            i_dft_end = i
            break
    dft_lines = lines[i_dft_beg:i_dft_end+1]
    
    in_structure = False
    for line in script_lines:
        if '&STRUCTURE' in line.upper(): in_structure = True

        if 'mol_subset_atoms' in line:
            continue #do not write this line
        elif 'molecule_axis' in line:
            settings_lines_frag.append(f'   molecule_axis = vector 1 0 0\n')
            continue
        else:
            found = False
            for key, val in settings_dict.items():
                if key == 'dft_settings_override' or val is None: continue
                if isinstance(val, list): val = ' '.join([str(x) for x in val])
                if key in line: 
                    settings_lines_frag.append(f'   {key} = {val}\n')
                    settings_dict.pop(key)
                    found = True
                    break
            if found: continue
            
        if in_structure and line.strip()=='/': 
            in_structure = False
            for key, val in settings_dict.items():
                if key == 'dft_settings_override' or not val: continue
                if isinstance(val, list): val = ' '.join([str(x) for x in val])
                settings_lines_frag.append(f'   {key} = {val}\n')

        settings_lines_frag.append(line)
    
    #add dft_settings_override
    settings_lines_frag += override_settings_adsorbed_fragment(program, dft_lines, settings_dict['dft_settings_override'].get('adsorption') if settings_dict['dft_settings_override'] else None)

    return settings_lines_frag


def setup_fragments_screening(RUN = False):
    
    with open(framgents_filename, "r") as f:
        fragments_dict = json.load(f)

    settings = Settings()
    
    main_dir = os.getcwd()
    for fragment_name, fragment_dict in fragments_dict["fragments"].items():
        
        os.chdir(main_dir)
        print(f'\n--------------------\nFragment {fragment_name}:')

        outdir = f'fragments/{fragment_name}'
        os.makedirs(outdir, exist_ok=True)

        shutil.copyfile(settings.jobscript, f'{outdir}/{jobscript_stdname}')

        frag_initial = generate_isolated_fragment(settings, fragment_dict)
        
        if not os.path.exists(FRAGMENTS_OUT_FILE_PATHS[settings.program].format(fragment_name)):
            print(f"Isolated {fragment_name} was not relaxed. Its initial structure will be used.")
            ase_custom.write_xyz_custom(f'{outdir}/{fragment_name}.xyz', frag_initial.mol_ase)
            fragment_filename = f'{outdir}/{fragment_name}.xyz'
        else:
            fragment_filename = FRAGMENTS_OUT_FILE_PATHS[settings.program].format(fragment_name)
    


        #only edit settings.in if not already present. This allows to fine-tune some parameters for specific fragments after -g, before -s
        answer = 'yes'
        if os.path.exists(f'{outdir}/settings.in'):
            print(f'{outdir}/settings.in already present.')
            while True:
                answer = input('Overwrite? ("y" = yes, "n" = no): ')
                if answer == 'yes' or answer == 'y' or answer == 'no' or answer == 'n': 
                    break
                else: print('Value not recognized. Try again.')
   
        if 'y' in answer:
            with open('settings.in', 'r') as f:
                settings_lines = f.readlines() 
            settings_dict = dict(jobscript=f'{jobscript_stdname} {settings.sbatch_command}',
                                slab_filename=os.path.abspath(settings.slab_filename),
                                molecule_filename=os.path.abspath(fragment_filename),
                                selected_atom_index=frag_initial.reference_atom_index,
                                x_rot_angles=fragment_dict.get("x_rot_angles", None),
                                y_rot_angles=fragment_dict.get("y_rot_angles", None),
                                z_rot_angles=fragment_dict.get("z_rot_angles", None),
                                vertical_angles=fragment_dict.get("vertical_angles", None),
                                screening_atom_distance=fragment_dict.get("screening_atom_distance", None),
                                screening_min_distance=fragment_dict.get("screening_min_distance", None),
                                relax_atom_distance=fragment_dict.get("relax_atom_distance", None),
                                relax_min_distance=fragment_dict.get("relax_min_distance", None),
                                fixed_indices_mol=frag_initial.constrained_indices,
                                fix_mol_xyz=fragment_dict.get("fix_mol_xyz", None),

                                dft_settings_override=fragment_dict.get("dft_settings_override", None)
                                )
            settings_lines = edit_fragment_settings(lines=settings_lines, settings_dict=settings_dict, program=settings.program)  
            with open(f'{outdir}/settings.in', 'w') as f:
                f.writelines(settings_lines)


        os.chdir(outdir)
        if RUN: launch_screening()
        os.chdir(main_dir)

    #print(json.dumps(fragments_dict, indent=3))


def final_relax_fragments(n_configs, threshold, BY_SITE):

    with open(framgents_filename, "r") as f:
        fragments_dict = json.load(f)

    main_dir = os.getcwd()
    for fragment_name in fragments_dict["fragments"]:

        os.chdir(main_dir)

        print(f'\n--------------------\nFragment {fragment_name}:')

        outdir = f'fragments/{fragment_name}'
        os.chdir(outdir)
        final_relax(n_configs=n_configs, threshold=threshold, BY_SITE=BY_SITE)
        os.chdir(main_dir)
        

def restart_jobs_fragments(calc_type : str):

    with open(framgents_filename, "r") as f:
        fragments_dict = json.load(f)

    main_dir = os.getcwd()
    for fragment_name in fragments_dict["fragments"]:

        os.chdir(main_dir)

        print(f'\n--------------------\nFragment {fragment_name}:')

        outdir = f'fragments/{fragment_name}'
        os.chdir(outdir)
        restart_jobs(calc_type)
        os.chdir(main_dir)


def get_diss_energies():

    with open(framgents_filename, "r") as f:
        fragments_dict = json.load(f)
    
    #get energy of most stable configuration for the whole molecule
    fragments_data = {"mol" : {}}

    datafile = pd.read_csv(labels_filename, index_col=0)
    settings = Settings(VERBOSE=False)
    results_mol = get_calculations_results(settings.program, 'RELAX', [0,0])
    energies_mol = []
    indices_mol = []
    for key, val in results_mol['energies'].items():
        energies_mol.append(val)
        indices_mol.append(key)
    i_min = indices_mol[energies_mol.index(min(energies_mol))]
    fragments_data["mol"]["energy"] = min(energies_mol)
    fragments_data["mol"]["site"] = datafile['site'][i_min]
    print("Mol energy collected.")


    #repeat for each fragment
    main_dir = os.getcwd()
    for fragment_name in fragments_dict["fragments"]:
        
        fragments_data.update({fragment_name : {} })

        
        outdir = f'fragments/{fragment_name}'
        os.chdir(outdir)
        datafile = pd.read_csv(labels_filename, index_col=0)
        results_mol = get_calculations_results(settings.program, 'RELAX', [0,0])
        energies_mol = []
        indices_mol = []
        for key, val in results_mol['energies'].items():
            energies_mol.append(val)
            indices_mol.append(key)
        i_min = indices_mol[energies_mol.index(min(energies_mol))]
        fragments_data[fragment_name]["energy"] = min(energies_mol)
        fragments_data[fragment_name]["site"] = datafile['site'][i_min]
        os.chdir(main_dir)

        print(f"{fragment_name} energy collected.")
    
    text = []

    for combination in fragments_dict["combinations"]:        
        
        products_names = combination[1]
        initial_fragment_name = combination[0]
        dissoc_products_total_ads_en = sum([ fragments_data[frag]["energy"] for frag in products_names])
        initial_fragment_energy = fragments_data[ initial_fragment_name ]["energy"]
        diss_en = dissoc_products_total_ads_en - initial_fragment_energy - settings.E_slab_mol[0] * (len(combination[1]) - 1)

        fragments_names = '{0}({1}) -> '.format(initial_fragment_name, fragments_data[initial_fragment_name]["site"]).format()    \
            +' + '.join([ '{0}({1})'.format(frag, fragments_data[frag]["site"]) for frag in products_names])

        text.append('{:70}{:+.3f}\n'.format(fragments_names, diss_en))

    with open('DISSOCIATION.txt', 'w') as f:
        f.write( '{0:70}{1}'.format("Fragments(most stable site)", "dissociation energy(eV)\n"))
        f.writelines( text )

    print(f"Dissociation results written to {'DISSOCIATION.txt'}")
