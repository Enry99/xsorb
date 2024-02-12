#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu 11 May 2023

@author: Enrico Pedretti

Generation of molecular fragments and dissociation energies calculation

"""

import json, os, shutil, copy
import numpy as np
from xsorbed.settings import Settings
from xsorbed.molecule import Molecule
from xsorbed.calculations import generate, final_relax
from xsorbed.io_utils import get_energies, restart_jobs
from xsorbed.common_definitions import *
from xsorbed.slab import Slab
from xsorbed.dftcode_specific import FRAGMENTS_OUT_FILE_PATHS, FRAGMENTS_IN_FILE_PATHS, SBATCH_POSTFIX_FRAGS, \
    override_settings_isolated_fragment, override_settings_adsorbed_fragment, Calculator
from xsorbed import ase_custom


TEST = True   #set to true for testing: prints sbatch command instead of actually launching jobs

#OK (code agnostic)
def generate_isolated_fragment(settings : Settings, fragment_dict : dict):
    """
    Generate isolated fragment from a molecule based on the provided settings and a dictionary for the fragment.

    Args:
       - settings (Settings): The settings object containing the required parameters for fragment generation.
       - fragment_dict (dict): A dictionary containing information about the fragment to be generated.

    Returns:
        list: A Molecule object representing the isolated fragment.
    """

        
    if "mol_subset_atoms" in fragment_dict and "break_bond_indices" in fragment_dict:
        raise ValueError("You can specify a fragment either by its atoms or by the broken bonds, not both.")
    
            
    axis = [x.strip("'") for x in fragment_dict['molecule_axis'].split()]
    if axis[0] == 'atoms':
        if len(axis[1:]) != 2:
            raise ValueError("Error: you must specify two integers as indices of the two atoms for the molecule axis.")
    elif axis[0] == 'vector':
        if len(axis[1:]) != 3:
            raise ValueError("Error: you must specify three numbers for the vector representing the molecule axis.")
    molecule_axis_atoms    = np.array(axis[1:], dtype=int).tolist() if axis[0] == 'atoms' else None
    axis_vector            = np.array(axis[1:], dtype=float).tolist() if axis[0] == 'vector' else None   
    
    if "mol_subset_atoms" in fragment_dict:
        mol = Molecule(molecule_filename=settings.molecule_filename,
                                atom_index=fragment_dict["selected_atom_index"],
                                molecule_axis_atoms=molecule_axis_atoms, 
                                axis_vector=axis_vector,
                                atoms_subset=fragment_dict["mol_subset_atoms"])
    elif "break_bond_indices" in fragment_dict:
        mol = Molecule(molecule_filename=settings.molecule_filename,
                                atom_index=fragment_dict["selected_atom_index"],
                                molecule_axis_atoms=molecule_axis_atoms, 
                                axis_vector=axis_vector,
                                break_bond_indices=fragment_dict["break_bond_indices"])
    else:
        raise RuntimeError("Error: fragment definition not recognized. You must specify either mol_subset_atoms or break_bond_indices.")
    

    return mol
    
#OK (code agnostic)
def generate_all_isolated_fragments(settings : Settings, fragments_dict : dict, VERBOSE : bool = True):
    
    structures = []
    for fragment_name, fragment_dict in fragments_dict["fragments"].items():

        if fragment_name == "mol":
            raise ValueError('Error: fragment name "mol" not allowed (reserved for whole molecule).')
        
        structures.append(generate_isolated_fragment(settings, fragment_dict))

    return structures

#OK (code agnostic) 
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
            override_settings_isolated_fragment(newsettings, molecules[0].Natoms, manual_dft_override)
        
        calc = Calculator(newsettings, fragment_name, molecule.mol_ase, outdir) 
        calc.write_input(molecule.mol_ase)

        frag_list.append(fragment_name)

    if INTERACTIVE: print('All input files written.') 

    return frag_list

#OK (code agnostic) 
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

#OK (code agnostic) 
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






def setup_fragments_screening(RUN = False, save_figs=False, saveas_format=None):
    
    with open(framgents_filename, "r") as f:
        fragments_dict = json.load(f)

    settings = Settings()
    
    main_dir = os.getcwd()
    
    for fragment_name, fragment_dict in fragments_dict["fragments"].items():
        
        os.chdir(main_dir)
        print(f'\n--------------------\nFragment {fragment_name}:')

        shutil.copyfile(settings.jobscript, f'fragments/{fragment_name}/{jobscript_stdname}')

        frag_initial = generate_isolated_fragment(settings, fragment_dict)
        
        if not os.path.exists(FRAGMENTS_OUT_FILE_PATHS[settings.program].format(fragment_name)):
            print(f"Isolated {fragment_name} was not relaxed. Its initial structure will be used.")
            mol = frag_initial
            #write to be able to use xsorb from within the fragment folder   
            ase_custom.write_xyz_custom(f'fragments/{fragment_name}/{fragment_name}.xyz', frag_initial.mol_ase) 

        else:
            mol = Molecule(molecule_filename=FRAGMENTS_OUT_FILE_PATHS[settings.program].format(fragment_name),
                        atom_index=frag_initial.reference_atom_index,
                        axis_vector=[1,0,0],  
                        fixed_indices_mol=frag_initial.constrained_indices, 
                        fix_mol_xyz=settings.fix_mol_xyz)
            

            #TODO:
            #mol_subset_atoms -> delete, since we are already providing the fragment as a input mol
            #molecule_axis -> set to vector 1 0 0, since the rotation is performed when generating the fragment
            #selected_atom_index -> reindex (just for completeness)
            #fixed_indices_mol -> reindex  (just for completeness)

            #TODO: check how to handle the dft parameters ordering.
    


        #only edit settings.in if not already present. This allows to fine-tune some parameters for specific fragments after -g, before -s
        answer = 'yes'
        if os.path.exists(f'fragments/{fragment_name}/settings.in'):
            print(f'fragments/{fragment_name}/settings.in already present.')
            while True:
                answer = input('Overwrite? ("y" = yes, "n" = no): ')
                if answer == 'yes' or answer == 'y' or answer == 'no' or answer == 'n': 
                    break
                else: print('Value not recognized. Try again.')

        if answer == 'yes' or answer == 'y': 
            #edit and put slab and molecule paths, and edit the jobscript changing the name to {jobscript_stdname}
    
            settings_lines_frag = copy.deepcopy(settings_lines)

            for i, line in enumerate(settings_lines_frag):
                if "molecule_filename" in line:
                    l = line.split('=')[:2]
                    l[1] = mol_frag_name.split('/')[-1]
                    line = '= '.join(l)+'\n'
                    settings_lines_frag[i] = line
    

            for flag in fragments_dict["fragments"][fragment_name]: 
                if flag ==  "mol_subset_atoms": continue
                found = False              
                for i, line in enumerate(settings_lines_frag):
                    
                    if flag in line:
                        l = line.split('=')[:2]
                        l[0] = '   '+flag.ljust(30)
                        if isinstance(fragments_dict["fragments"][fragment_name][flag], list):
                            l[1] = ' '.join(str(x) for x in fragments_dict["fragments"][fragment_name][flag])                        
                        else:
                            l[1] = str(fragments_dict["fragments"][fragment_name][flag])
                        line = '= '.join(l)+'\n'
                        settings_lines_frag[i] = line
                        found = True
                        break
                if not found:
                    if isinstance(fragments_dict["fragments"][fragment_name][flag], list):
                        l = ' '.join(str(x) for x in fragments_dict["fragments"][fragment_name][flag])                        
                    else:
                        l = str(fragments_dict["fragments"][fragment_name][flag])
                    settings_lines_frag.insert(STRUCTURE_line_i, '   '+flag.ljust(30)+' = '+l+'\n')


            with open('fragments/{0}'.format(fragment_name)+'/settings.in', 'w') as f:
                f.writelines(settings_lines_frag)


        
        j_dir = 'fragments/'+fragment_name
        if not os.path.exists(j_dir + '/'+jobscript_filename) : shutil.copyfile(settings.jobscript, j_dir + '/'+jobscript_filename)

        os.chdir(j_dir)
        generate(RUN = RUN, etot_forc_conv = etot_forc_conv, SAVEFIG=save_figs, saveas_format=saveas_format)
        os.chdir(main_dir)

    #print(json.dumps(fragments_dict, indent=3))


def final_relax_fragments(n_configs, threshold, BY_SITE):

    with open(framgents_filename, "r") as f:
        fragments_dict = json.load(f)

    main_dir = os.getcwd()
    for fragment_name in fragments_dict["fragments"]:

        os.chdir(main_dir)

        print('\n--------------------\nFragment {}:'.format(fragment_name))

        j_dir = 'fragments/'+fragment_name
        os.chdir(j_dir)
        final_relax(n_configs=n_configs, threshold=threshold, BY_SITE=BY_SITE)
        os.chdir(main_dir)
        

def restart_jobs_fragments(which):

    with open(framgents_filename, "r") as f:
        fragments_dict = json.load(f)

    main_dir = os.getcwd()
    for fragment_name in fragments_dict["fragments"]:

        os.chdir(main_dir)

        print('\n--------------------\nFragment {}:'.format(fragment_name))

        j_dir = 'fragments/'+fragment_name
        os.chdir(j_dir)
        restart_jobs(which='screening' if which == 's' else 'relax')
        os.chdir(main_dir)


def get_diss_energies():

    from natsort import natsorted
    import glob

    with open(framgents_filename, "r") as f:
        fragments_dict = json.load(f)
    
    
    slab_energy = Settings().E_slab_mol[0] #* rydbergtoev  CHANGE!!!


    fragments_data = {"mol" : {}}
    #get total energy of molecule, each with the label and position
    sites = [] 
    config_labels = []
    with open('site_labels.csv', 'r') as f:
        file = f.readlines()
        for line in file:
            if 'site' in line: continue
            sites.append(int(line.split(',')[0]))
            config_labels.append(line.split(',')[4])
    energies = get_energies(pwo_prefix='relax')
    i_min = energies.index(min(energies))
    files = natsorted(glob.glob( 'relax_*.pwo' ))
    name_min = files[i_min]
    label_min = int(name_min.split('_')[-1].split('.')[0])
    fragments_data["mol"]["energy"] = min(energies)
    fragments_data["mol"]["site"] = config_labels[label_min]

    #repeat for each fragment
    for fragment_name in fragments_dict["fragments"]:
        
        fragments_data.update({fragment_name : {} })

        #get total energy of each fragment, each with the label and position
        sites = [] 
        config_labels = []
        with open('fragments/{}/site_labels.csv'.format(fragment_name), 'r') as f:
            file = f.readlines()

            for line in file:
                if 'site' in line: continue
                sites.append(int(line.split(',')[0]))
                config_labels.append(line.split(',')[4])

        energies = get_energies(pwo_prefix='fragments/{}/relax'.format(fragment_name))

        i_min = energies.index(min(energies))
        
        files = natsorted(glob.glob( 'fragments/{}/relax_*.pwo'.format(fragment_name) ))

        name_min = files[i_min]
        label_min = int(name_min.split('_')[-1].split('.')[0])

        fragments_data[fragment_name]["energy"] = min(energies)
        fragments_data[fragment_name]["site"] = config_labels[label_min]
        

    
    text = []


    for combination in fragments_dict["combinations"]:        
        
        products_names = combination[1]
        initial_fragment_name = combination[0]
        dissoc_products_toten = sum([ fragments_data[frag]["energy"] for frag in products_names])
        initial_fragment_energy = fragments_data[ initial_fragment_name ]["energy"]
        diss_en = dissoc_products_toten - initial_fragment_energy - slab_energy * (len(products_names) - 1)

        fragments_names = '{0}({1}) -> '.format(initial_fragment_name, fragments_data[initial_fragment_name]["site"]).format()    \
            +' + '.join([ '{0}({1})'.format(frag, fragments_data[frag]["site"]) for frag in products_names])

        text.append('{:70}{:+.3f}\n'.format(fragments_names, diss_en))

    with open('DISSOCIATION.txt', 'w') as f:
        f.write( '{0:70}{1}'.format("Fragments(most stable site)", "dissociation energy(eV)\n"))
        f.writelines( text )
