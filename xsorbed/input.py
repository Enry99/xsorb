#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue 28 Feb 2023
@author: Enrico Pedretti

Small module with function definitions to read settings from input file
and return them as two dictionaries (script settings and Espresso settings)

"""

import sys

#NOTE 1: The blocks CELL_PARAMETERS ATOMIC_POSITIONS ATOMIC_SPECIES must NOT be included in input file, as they are read from the input structures
#
#NOTE 3(possible TODO): This code does not yet support the following Espresso blocks:
#OCCUPATIONS, CONSTRAINTS, ATOMIC_VELOCITIES, ATOMIC_FORCES, ADDITIONAL_K_POINTS, SOLVENTS, HUBBARD

#functions
def _is_number(s : str):
    try:
        float(s)
        return True
    except ValueError:
        return False


def _read_card(lines : list[str], CONVERT = False):

    block_dict =  {}

    for line in lines:

        (key, val) = line.split('=')

        key = key.strip()
        val = val.strip()
        val = val.strip("'")

        if(CONVERT):
            if (_is_number(val)):
                if val.isnumeric(): val = int(val)
                else: val = float(val)      

        block_dict[key] = val

    return  block_dict


def read_input_file(filename: str):

    script_settings_dict   = {}    
    espresso_settings_dict = {}


    with open(filename) as file:
        
        lines = []
        kpoints = []  #syntax: ['gamma'/'automatic', [kx,ky,kz], [koffx,koffy,koffz]]
        atomic_species = []  #syntax: [[Fe, 1.000, Fe.pbe-n-rrkjus_psl.1.0.0.UPF], ...]
        last_dump = [] #other lines not in the two main blocks
        in_card = False
        ATOMIC_SPECIES = False
        KPOINTS = False
        espresso_block = False
        settings_block = False

        for line in file:
            
            if line.isspace() or line.strip()[0] == '!' or line.strip()[0] == '#':
                if(ATOMIC_SPECIES): ATOMIC_SPECIES = False
                continue #skip empty / comment lines
            if(ATOMIC_SPECIES):
                if ('ATOMIC_POSITIONS' in line or 
                'K_POINTS' in line or 
                'ADDITIONAL_K_POINTS' in line or 
                'CELL_PARAMETERS' in line or 
                'CONSTRAINTS' in line or 
                'OCCUPATIONS' in line or 
                'ATOMIC_VELOCITIES' in line or 
                'ATOMIC_FORCES' in line or 
                'SOLVENTS' in line or 
                'HUBBARD' in line):
                    ATOMIC_SPECIES = False

            
            if('#' in line):
                line = line.split('#')[0] #deal with comments in the lines
            if('!' in line):
                line = line.split('!')[0]
            
            if "@ESPRESSO" in line:
                espresso_block = True
                continue
            if "@SETTINGS" in line:
                settings_block = True
                continue

            if "@/ESPRESSO" in line:
                espresso_block = False
                continue
            if "@/SETTINGS" in line:
                settings_block = False
                continue


            if(settings_block):
                if line.strip()[0] == '&':
                    block_name = line.split('&')[1].strip()
                    in_card = True
                    continue

                if line.strip()=='/':
                    script_settings_dict.update({block_name : _read_card(lines)})
                    in_card = False
                    lines.clear()
                    continue

                if in_card: lines.append(line)  

            elif(espresso_block):

                if line.strip()[0] == '&':
                    block_name = line.split('&')[1].strip()
                    in_card = True
                    continue

                if line.strip()=='/':
                    espresso_settings_dict.update({block_name : _read_card(lines, CONVERT=True)})
                    in_card = False
                    lines.clear()
                    continue

                if in_card: 
                    lines.append(line)
                    continue



                if('ATOMIC_SPECIES' in line or 'atomic_species' in line):
                    ATOMIC_SPECIES = True
                    continue
                if(ATOMIC_SPECIES):
                    atomic_species.append(line.split())
                    continue

                if('K_POINTS' in line or 'kpoints' in line):
                    KPOINTS = True
                    gamma_or_auto = line.split()[1].strip().lower()
                    kpoints.append(gamma_or_auto)
                    if 'gamma' in gamma_or_auto:
                        KPOINTS = False
                    continue
                if(KPOINTS):
                    kpoints.append(line.split()[:3])
                    kpoints.append(line.split()[3:])
                    KPOINTS = False
                    continue

                #collects everything that did not fall in any of the previous category. This is appended as-is at the end of the pwi.
                last_dump.append(line)         
    

    if(not script_settings_dict):
        print('Script settings not read correctly. Quitting.')
        sys.exit(1)
    if(not espresso_settings_dict):
        print('Espresso settings not read correctly. Quitting.')
        sys.exit(1)
    if(not kpoints):
        print('Kpoints not read correctly. Quitting.')
        sys.exit(1)
    if(not atomic_species):
        print('Atomic_species not read correctly. Quitting.')
        sys.exit(1)

    return script_settings_dict, espresso_settings_dict, atomic_species, kpoints, last_dump
