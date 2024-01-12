#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu 11 May 2023

@author: Enrico Pedretti

Generation of molecular fragments and dissociation energies calculation

"""
#TODO: fragments, se atomo isolato aggiungere nosym=.true.

import json, os, shutil, copy, sys
from settings import Settings
from molecule import Molecule
from espresso_mod import Espresso_mod
from calculations import generate, final_relax
from io_utils import get_energies, restart_jobs
from filenames import *
from slab import Slab


TEST = False   #set to true for testing: prints sbatch command instead of actually launching jobs


def isolated_fragments(RUN = False):

    #NOTE!!! By default the program sets starting magnetizations to 1 and mixing_beta to 0.1 for all fragments to facilitate convergence. Might be changed in the future.


    with open("fragments.json", "r") as f:
        fragments_dict = json.load(f)

    settings = Settings()

    mol_whole = Molecule(settings.molecule_filename)

    if not os.path.exists('fragments'): os.mkdir('fragments')
    for fragment_name in fragments_dict["fragments"]:

        print("Generating fragment {}...".format(fragment_name))

        if fragment_name == "mol":
            print('Error: fragment name "mol" not allowed (reserved for whole molecule). Quitting.')
            sys.exit(1)

        if not os.path.exists('fragments/'+fragment_name): os.mkdir('fragments/'+fragment_name)


        if os.path.exists('fragments/'+fragment_name+'/'+fragment_name+'.pwo'): 
            print('fragments/'+fragment_name+'/'+fragment_name+'.pwo already present, possibly from a running calculation. You can decide to re-calculate it or skip it.')        
            while True:
                answer = input('Re-calculate? ("y" = yes, "n" = no): ')
                if answer == 'yes' or answer == 'y' or answer == 'no' or answer == 'n': 
                    break
                else: print('Value not recognized. Try again.')
            if 'n' in answer: continue

        answer = 'yes'
        if os.path.exists('fragments/'+fragment_name+'/'+fragment_name+'.pwi'):
            print(fragment_name+'.pwi already present.')
            while True:
                answer = input('Overwrite? ("y" = yes, "n" = no): ')
                if answer == 'yes' or answer == 'y' or answer == 'no' or answer == 'n': 
                    break
                else: print('Value not recognized. Try again.')

        if answer == 'yes' or answer == 'y':     
            #create complementary
            if "mol_subset_atoms" in fragments_dict["fragments"][fragment_name] and "mol_subset_atoms_compl" in fragments_dict["fragments"][fragment_name]:
                print("You can specify a fragment either by its atoms or by the complementary atoms of the whole molecule, not both. Quitting.")
                sys.exit(1)
            
            if "mol_subset_atoms_compl" in fragments_dict["fragments"][fragment_name]:
                mol_subset_atoms = [i for i in range(mol_whole.natoms) if i not in fragments_dict["fragments"][fragment_name]["mol_subset_atoms"]]
            else:
                mol_subset_atoms = fragments_dict["fragments"][fragment_name]["mol_subset_atoms"]
            
            mol = Molecule(settings.molecule_filename, atoms_subset=mol_subset_atoms)      

            fragment_species = set([atom.symbol for atom in mol.mol_ase])            
            absent_species = [species for species in settings.pseudopotentials if species not in fragment_species]

            pseudos_subset = copy.deepcopy(settings.pseudopotentials)
            for species in absent_species:
                pseudos_subset.pop(species)


            absent_indices = []
            for i, pot in enumerate(settings.pseudopotentials):
                if pot in absent_species: absent_indices.append(i)

            starting_mag_subset = []
            for i, mag in enumerate(settings.starting_mag):
                if i not in absent_indices:
                    #starting_mag_subset.append(mag)
                    starting_mag_subset.append(1.0) #set all to 1, can be better for fragments
            
            delta_n = []
            missing_counter = 0
            for j, species in enumerate(settings.pseudopotentials):
                if species in absent_species:
                    missing_counter += 1
                delta_n.append(missing_counter)


            flags_i_subset      = []
            settings_flags_i = copy.deepcopy(settings.flags_i)
            for flag in settings_flags_i:
                i_flag = int(flag[0].split('(')[1].split(')')[0])
                if i_flag-1 not in absent_indices:
                    #print(i_flag)
                    flag[0] = flag[0].replace('({0})'.format(i_flag), '({0})'.format(i_flag - delta_n[i_flag-1]))
                    flags_i_subset.append(flag)

            settings.espresso_settings_dict['ELECTRONS'].update({'mixing_beta' : 0.1})
            if(len(mol.mol_ase) == 1): settings.espresso_settings_dict['SYSTEM'].update({'nosym' : True}) #isolated atoms

            calc = Espresso_mod(pseudopotentials=pseudos_subset, 
                    input_data=settings.espresso_settings_dict,
                    filename='fragments/'+fragment_name+'/'+fragment_name+'.pwi',
                    kpts = None,   #TODO: handle case where gamma is not right, here we assume that the input cell for the molecule was big enough
                    koffset= None)
            calc.set_system_flags(starting_mag_subset, flags_i_subset)
            calc.write_input(mol.mol_ase)


        if(RUN):
            sbatch_command = settings.sbatch_command 

            main_dir = os.getcwd()

            input_file = fragment_name+'.pwi'
            output_file = input_file.replace("pwi", "pwo")


            j_dir = 'fragments/'+fragment_name
            shutil.copyfile(settings.jobscript, j_dir + '/'+jobscript_filename)

            os.chdir(j_dir) #####################
            with open(jobscript_filename, 'r') as f:
                lines = f.readlines()

                for i, line in enumerate(lines):
                    if "job-name" in line:
                        lines[i] = line.split('=')[0] + '="' + fragment_name + '"\n'
                        break
        
            with open(jobscript_filename, 'w') as f:
                f.writelines( lines )


            if(TEST): print(sbatch_command+" " + jobscript_filename + ' ' +input_file + ' '+output_file)
            else: os.system(sbatch_command+" " + jobscript_filename + ' ' +input_file + ' '+output_file)  #launchs the jobscript in j_dir from j_dir
            os.chdir(main_dir) ####################


def setup_fragments_screening(RUN = False, etot_forc_conv = hybrid_screening_thresholds, save_figs=False, saveas_format=None):
    with open("fragments.json", "r") as f:
        fragments_dict = json.load(f)

    settings = Settings()
    slab = Slab(settings.slab_filename)  
 
    settings_lines = settings.text

    for i, line in enumerate(settings_lines):
        if '&STRUCTURE' in line:
            STRUCTURE_line_i = i+1
            break


    for fragment_name in fragments_dict["fragments"]:
        if not os.path.exists('fragments/{0}/{0}.pwi'.format(fragment_name)) and not os.path.exists('fragments/{0}/{0}.pwo'.format(fragment_name)):
            print("ERROR! fragments/{0}/{0} (.pwi or .pwo) not found. You need to generate/relax all the fragments first.".format(fragment_name))
            sys.exit(1)

    main_dir = os.getcwd()

    for fragment_name in fragments_dict["fragments"]:

        os.chdir(main_dir)

        print('\n--------------------\nFragment {}:'.format(fragment_name))

        mol_frag_name = "fragments/{0}/{0}.{1}".format(fragment_name, 'pwo' if os.path.exists('fragments/{0}/{0}.pwo'.format(fragment_name)) else 'pwi')
        mol = Molecule(mol_frag_name)
    
        if not os.path.exists('fragments/{0}/{1}'.format(fragment_name, settings.slab_filename)): shutil.copyfile(settings.slab_filename, 'fragments/{0}/{1}'.format(fragment_name, settings.slab_filename))

        shutil.copyfile(settings.jobscript, 'fragments/{0}/{1}'.format(fragment_name, jobscript_filename))

        #only edit settings.in if not already present. This allows to fine-tune some parameters for specific fragments after -g, before -s

        answer = 'yes'
        if os.path.exists('fragments/{0}/settings.in'.format(fragment_name)):
            print('fragments/{0}/settings.in already present.'.format(fragment_name))
            while True:
                answer = input('Overwrite? ("y" = yes, "n" = no): ')
                if answer == 'yes' or answer == 'y' or answer == 'no' or answer == 'n': 
                    break
                else: print('Value not recognized. Try again.')

        if answer == 'yes' or answer == 'y': 
    
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


            slab_species = set(slab.slab_ase.get_chemical_symbols())
            fragment_species = set([atom.symbol for atom in mol.mol_ase])            
            absent_species = [species for species in settings.pseudopotentials if (species not in fragment_species and species not in slab_species)]
            absent_indices = []
            for i, pot in enumerate(settings.pseudopotentials):
                if pot in absent_species: absent_indices.append(i+1)


            #remove undesired lines
            pseudos = list(settings.pseudopotentials.keys())
            for species in absent_species:
                for i, line in enumerate(settings_lines_frag):
                    if settings.pseudopotentials[species] in line:
                        del settings_lines_frag[i]
                    if '({0})'.format(pseudos.index(species)+1) in line:
                        del settings_lines_frag[i]

            #reindex flags of the format (i)
            missing_counter = 0
            for j, species in enumerate(pseudos):
                if species in absent_species:
                    missing_counter += 1
                    continue
                for i, line in enumerate(settings_lines_frag):           
                    if '({0})'.format(j+1) in line:
                        settings_lines_frag[i] = line.replace('({0})'.format(j+1), '({0})'.format(j+1 - missing_counter))

            with open('fragments/{0}'.format(fragment_name)+'/settings.in', 'w') as f:
                f.writelines(settings_lines_frag)


        
        j_dir = 'fragments/'+fragment_name
        if not os.path.exists(j_dir + '/'+jobscript_filename) : shutil.copyfile(settings.jobscript, j_dir + '/'+jobscript_filename)

        os.chdir(j_dir)
        generate(RUN = RUN, etot_forc_conv = etot_forc_conv, SAVEFIG=save_figs, saveas_format=saveas_format)
        os.chdir(main_dir)

    #print(json.dumps(fragments_dict, indent=3))


def final_relax_fragments(n_configs, threshold, BY_SITE):

    with open("fragments.json", "r") as f:
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

    with open("fragments.json", "r") as f:
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

    with open("fragments.json", "r") as f:
        fragments_dict = json.load(f)
    
    
    slab_energy = Settings().E_slab_mol[0] * rydbergtoev


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
