#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu 11 May 2023

@author: Enrico Pedretti

Generation of molecular fragments and dissociation energies calculation

"""


import json, os, shutil, copy, sys
from settings import Settings
from molecule import Molecule
from espresso_mod import Espresso_mod
from calculations import generate, final_relax
from io_utils import get_energies, restart_jobs
from filenames import *
from slab import Slab


TEST = False


def isolated_fragments(RUN = False):
    with open("fragments.json", "r") as f:
        fragments_dict = json.load(f)

    settings = Settings()

    mol_whole = Molecule(settings.molecule_filename)

    if not os.path.exists('fragments'): os.mkdir('fragments')
    for fragment_name in fragments_dict["fragments"]:

        if not os.path.exists('fragments/'+fragment_name): os.mkdir('fragments/'+fragment_name)
        
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
                starting_mag_subset.append(mag)
        
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

            if(os.path.isfile(output_file)): 
                print(output_file+' already present, possibly from a running calculation. It will be skipped.')
                continue

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


    if not os.path.exists('fragments'): os.mkdir('fragments')

    #TODO: check if ALL fragments exist in folders, otherwise ERROR! you need to generate/relax the fragments first

    main_dir = os.getcwd()

    for fragment_name in fragments_dict["fragments"]:

        os.chdir(main_dir)

        print('\n--------------------\nFragment {}:'.format(fragment_name))

        mol_frag_name = "'{0}.{1}'".format(fragment_name, 'pwo' if os.path.exists('fragments/{0}/{0}.pwo'.format(fragment_name)) else 'pwi')
        mol = Molecule(mol_frag_name)
    
        if not os.path.exists('fragments/{0}/{1}'.format(fragment_name, settings.slab_filename)): shutil.copyfile(settings.slab_filename, 'fragments/{0}/{1}'.format(fragment_name, settings.slab_filename))

        shutil.copyfile(settings.jobscript, 'fragments/{0}/{1}'.format(fragment_name, jobscript_filename))

        #only edit settings.in if not already present. This allows to fine-tune some parameters for specific fragments after -g, before -s
        if not os.path.exists('fragments/{0}/settings.in'.format(fragment_name)): 
    
            settings_lines_frag = copy.deepcopy(settings_lines)

            for i, line in enumerate(settings_lines_frag):
                if "molecule_filename" in line:
                    l = line.split('=')[:2]
                    l[1] = mol_frag_name
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
    
    
    whole_mol_energy = get_energies(pwo_prefix='relax')

    fragments_data = {}
    for fragment_name in fragments_dict["fragments"]:
        
        fragments_data.update({fragment_name : {} })

        #get total energy of each fragment, each with the label and position
        site_labels = []
        with open('fragments/{}/site_labels.csv'.format(fragment_name), 'r') as f:
            file = f.readlines()

            for line in file:
                if 'site' in line: continue
                site_labels.append(int(line.split(',')[4].split(' ')[0]))

        energies = get_energies(pwo_prefix='fragments/{}/relax'.format(fragment_name))

        i_min = energies.index(min(energies))
        
        files = natsorted(glob.glob( 'fragments/{}/relax'.format(fragment_name) ))

        name_min = files[i_min]
        label_min = int(name_min.split('_')[-1])

        fragments_data[fragment_name]["energy"] = min(energies)
        fragments_data[fragment_name]["site"] = site_labels[label_min]
        

    
    text = []


    for combination in fragments_dict["combinations"]:        
        
        fragments_toten = sum([ fragments_data[frag]["energy"] for frag in combination])
        diss_en = fragments_toten - whole_mol_energy

        fragments_names = ' + '.join([ '{0}({1})'.format(frag, fragments_data[frag]["site"]) for frag in combination])

        text.append('{:50}{:.3f}\n'.format(fragments_names, diss_en))

    with open('dissociation.txt', 'w') as f:
        f.write( '{0:50}{1}'.format("Fragments(most stable site)", "dissociation energy(eV)\n"))
        f.writelines( text )
