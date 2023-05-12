#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu 11 May 2023

@author: Enrico Pedretti

Helper to handle generation of molecular fragments

"""

import json, os
from settings import Settings
from molecule import Molecule
from slab import Slab


def read_fragments():
    with open("fragments.json", "r") as f:
        fragments_dict = json.load(f)

    settings = Settings()
    slab = Slab(settings.slab_filename)
    mol = Molecule(settings.molecule_filename)
 
    settings_lines = settings.text

    for i, line in enumerate(settings_lines):
        if '&STRUCTURE' in line:
            STRUCTURE_line_i = i+1
            break


    if not os.path.exists('fragments'): os.mkdir('fragments')

    for fragment_name in fragments_dict.keys():
        if not os.path.exists('fragments/'+fragment_name): os.mkdir('fragments/'+fragment_name)

        names_AB = list(fragments_dict[fragment_name].keys())

        fragments_dict[fragment_name][names_AB[1]]["mol_subset_atoms"] = [i for i in range(mol.natoms) if i not in fragments_dict[fragment_name][names_AB[0]]["mol_subset_atoms"]]

        for name in names_AB:
            if not os.path.exists('fragments/{0}/{1}'.format(fragment_name, name)): os.mkdir('fragments/{0}/{1}'.format(fragment_name, name))

            settings_lines_AB = settings_lines.copy()


            for flag in fragments_dict[fragment_name][name]:  
                found = False              
                for i, line in enumerate(settings_lines_AB):
                    
                    if flag in line:
                        l = line.split('=')[:2]
                        l[0] = '   '+flag.ljust(30)
                        if isinstance(fragments_dict[fragment_name][name][flag], list):
                            l[1] = ' '.join(str(x) for x in fragments_dict[fragment_name][name][flag])                        
                        else:
                            l[1] = str(fragments_dict[fragment_name][name][flag])
                        line = '= '.join(l)+'\n'
                        settings_lines_AB[i] = line
                        found = True
                        break
                if not found:
                    if isinstance(fragments_dict[fragment_name][name][flag], list):
                        l = ' '.join(str(x) for x in fragments_dict[fragment_name][name][flag])                        
                    else:
                        l = str(fragments_dict[fragment_name][name][flag])
                    settings_lines_AB.insert(STRUCTURE_line_i, '   '+flag.ljust(30)+' = '+l+'\n')


            slab_species = set(slab.slab_ase.get_chemical_symbols())
            fragment_species = set([atom.symbol for atom in mol.mol_ase if atom.index in fragments_dict[fragment_name][name]["mol_subset_atoms"]])            
            absent_species = [species for species in settings.pseudopotentials if (species not in fragment_species and species not in slab_species)]
            absent_indices = []
            for i, pot in enumerate(settings.pseudopotentials):
                if pot in absent_species: absent_indices.append(i+1)


            #remove undesired lines
            pseudos = list(settings.pseudopotentials.keys())
            for species in absent_species:
                for i, line in enumerate(settings_lines_AB):
                    if settings.pseudopotentials[species] in line:
                        del settings_lines_AB[i]
                    if '({0})'.format(pseudos.index(species)+1) in line:
                        del settings_lines_AB[i]

            #reindex flags of the format (i)
            j = 0
            for i, line in enumerate(settings_lines_AB):
                for index in absent_indices:
                    j-=1
                    if '({0})'.format(i) in line:
                        settings_lines_AB[i].replace('({0})'.format(index), '({0})'.format(index-j))

            with open('fragments/{0}/{1}'.format(fragment_name, name)+'/settings.in', 'w') as f:
                f.writelines(settings_lines_AB)

    #Settings().    

    #print(json.dumps(fragments_dict, indent=3))

