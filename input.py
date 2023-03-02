#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue 28 Feb 2023
@author: Enrico Pedretti

Small module with function definitions to read settings from input file
and return them as two dictionaries (script settings and Espresso settings)

"""


#NOTE 1: The blocks CELL_PARAMETERS ATOMIC_POSITIONS ATOMIC_SPECIES must NOT be included in input file, as they are read from the input structures
#NOTE 2: The kpoints, kpoints offsets and coordinates to be fixed must be set in the first section "@SETTINGS" of the input file, not in "@ESPRESSO"
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


def _read_block(lines : list[str], CONVERT = False):

    block_dict =  {}

    for line in lines:

        (key, val) = line.split('=')

        key = key.strip()
        val = val.strip()
        val = val.strip("'")

        if(CONVERT):
            if (_is_number(val)): #don't do the conversion here, but in the main module
                if val.isnumeric(): val = int(val)
                else: val = float(val)      

        block_dict[key] = val

    return  block_dict


def read_input_file(filename: str):

    script_settings_dict   = {}    
    espresso_settings_dict = {}


    with open(filename) as file:
        
        lines = []
        last_dump = [] #other lines not in the two main blocks
        in_block = False
        espresso_block = False
        settings_block = False

        for line in file:
            
            if line.isspace() or line.strip()[0] == '!' or line.strip()[0] == '#': continue #skip empty / comment lines
            
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
                    in_block = True
                    continue

                if line.strip()=='/':
                    script_settings_dict.update({block_name : _read_block(lines)})
                    in_block = False
                    lines.clear()
                    continue

                if in_block: lines.append(line)  

            elif(espresso_block):

                if line.strip()[0] == '&':
                    block_name = line.split('&')[1].strip()
                    in_block = True
                    continue

                if line.strip()=='/':
                    espresso_settings_dict.update({block_name : _read_block(lines, CONVERT=True)})
                    in_block = False
                    lines.clear()
                    continue

                if in_block: lines.append(line)
                else: last_dump.append(line)
    

    if(not script_settings_dict):
        raise RuntimeError('Script settings not read correctly. Quitting.')
    if(not espresso_settings_dict):
        raise RuntimeError('Espresso settings not read correctly. Quitting.')

    return script_settings_dict, espresso_settings_dict
