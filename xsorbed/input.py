#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue 28 Feb 2023
@author: Enrico Pedretti

Small module with function definitions to read settings from input file
and return them as two dictionaries (script settings and dft program settings)

"""

from dftcode_specific import get_dftprogram_settings, SUPPORTED_PROGRAMS


def _is_number(s : str):
    try:
        float(s)
        return True
    except ValueError:
        return False


def _read_card(lines : list, CONVERT : bool = False):
    '''
    Read a card with Espresso-style formatting, returning a dictionary with the entries

    Args:
    - lines: list of strings containing the card content
    - CONVERT: if True, try to convert string values to numbers
    '''

    card_dict =  {}

    for line in lines:

        (key, val) = line.split('=')

        key = key.strip()
        val = val.strip()
        val = val.strip("'")
        val = val.strip('"')

        key = key.lower()
        #val = val.lower()

        if(CONVERT):
            if (_is_number(val)):
                if val.isnumeric(): val = int(val)
                else: val = float(val)
            else:
                if 'true' in val.lower(): val=True
                elif 'false' in val.lower(): val=False   

        card_dict[key] = val

    return  card_dict


def read_input_file(filename: str):
    '''
    Parse the settings.in file, returning a dictionary with the STRUCTURE settings, 
    and a dictionary with the settings for the DFT code

    Args:
    - filename: path of the settings.in file

    Returns:
    - script_settings_dict: dictionary with the settings for adsorption structures generation. In the dictionary
    the values are here ALL inserted as strings, and need to be converted later to the appropriate format
    - dftprogram_settings_dict: dictionary with the specific settings for the dft program
    - PROGRAM: name of the DFT program
    '''

    script_settings_dict   = {}    
    dftprogram_settings_dict = {}

    with open(filename) as file:
        
        dftprogram_lines = []
        lines = []
        in_card = False
        settings_block = False
        dftprogram_block = False

        for line in file:

            #skip empty / comment lines
            if line.isspace() or line.strip()[0] == '!' or line.strip()[0] == '#':
                continue 

            #deal with comments in the lines
            if('#' in line): line = line.split('#')[0] 
            if('!' in line): line = line.split('!')[0]            
            
                
            #check which block we are in and sets to True/False the corresponding variable
            for prog in SUPPORTED_PROGRAMS: 
                if f"@{prog}" in line.upper():
                    PROGRAM = prog
                    dftprogram_block = True
                    continue
            
            if dftprogram_block and f"@/{PROGRAM}" in line.upper():
                dftprogram_block = False
                continue

            if "@SETTINGS" in line.upper():
                settings_block = True
                continue

            if "@/SETTINGS" in line.upper():
                settings_block = False
                continue


            #read the block content
            if(settings_block):
                
                #begin card
                if line.strip()[0] == '&':
                    block_name = line.split('&')[1].strip().upper()
                    in_card = True
                    continue

                #end card, add collected content to dict
                if line.strip()=='/':
                    script_settings_dict.update({block_name : _read_card(lines)})
                    in_card = False
                    lines.clear()
                    continue

                if in_card: lines.append(line)  

            elif(dftprogram_block):
                dftprogram_lines.append(line)
                continue
 
    
    dftprogram_settings_dict = get_dftprogram_settings(PROGRAM, dftprogram_lines)

    if(not script_settings_dict):
        raise RuntimeError('Script settings not read correctly.')
    if(not dftprogram_settings_dict):
        raise RuntimeError(f'{PROGRAM} settings not read correctly.')

    return script_settings_dict, dftprogram_settings_dict, PROGRAM
