#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Enrico Pedretti

'''
Module to read input settings for different DFT programs and store them in a dictionary.
'''

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional
import sys

from dacite import from_dict
from ase.io.espresso import read_fortran_namelist

from xsorb.dft_codes.definitions import HYBRID_SCREENING_THRESHOLDS
from xsorb.io import settings

# dataclasses to store DFT program settings
@dataclass
class EspressoSettings:
    '''
    Dataclass to store Espresso input settings.
    '''
    pwi_path: str
    pwi_path_screening : Optional[str]
    etot_conv_thr_screening: float = HYBRID_SCREENING_THRESHOLDS['espresso'][0]
    forc_conv_thr_screening: float = HYBRID_SCREENING_THRESHOLDS['espresso'][1]

    # non-initialized attributes (filled in __post_init__)
    settings_dict = None

    def __post_init__(self):

        #NOTE 1: The blocks CELL_PARAMETERS ATOMIC_POSITIONS ATOMIC_SPECIES must NOT be included
        # in input file, as they are read from the input structures
        #NOTE 2: This code does not yet support the following Espresso blocks:
        #OCCUPATIONS, CONSTRAINTS, ATOMIC_VELOCITIES, ATOMIC_FORCES, ADDITIONAL_K_POINTS, SOLVENTS

        # parse namelist section and extract remaining lines
        with open(self.pwi_path, 'r') as file:
            settings_dict, card_lines = read_fortran_namelist(file)
            self.settings_dict = dict(settings_dict)


        #parse ATOMIC_SPECIES, K_POINTS and HUBBARD
        atomic_species_index, k_points_index, hubbard_index = None, None, None
        for i, line in enumerate(card_lines):
            if 'ATOMIC_SPECIES' in line.upper():
                atomic_species_index = i
            if 'K_POINTS' in line.upper():
                k_points_index = i
            if 'HUBBARD' in line.upper():
                hubbard_index = i
        if atomic_species_index is None:
            raise ValueError('ATOMIC_SPECIES card not found in input file.')
        if k_points_index is None:
            raise ValueError('K_POINTS card not found in input file.')

        def _end_of_card(card, line):

            cards_list = ['ATOMIC_SPECIES',
                        'ATOMIC_POSITIONS',
                        'K_POINTS',
                        'ADDITIONAL_K_POINTS',
                        'CELL_PARAMETERS',
                        'CONSTRAINTS',
                        'OCCUPATIONS',
                        'ATOMIC_VELOCITIES',
                        'ATOMIC_FORCES',
                        'SOLVENTS',
                        'HUBBARD']

            for other_card in cards_list:
                if other_card == card.upper(): continue
                if other_card in line.upper(): return True

            return False

        #ATOMIC_SPECIES
        self.settings_dict['pseudopotentials'] = {}
        i = atomic_species_index+1
        while i < len(card_lines):
            line = card_lines[i]

            if _end_of_card(card='ATOMIC_SPECIES', line=line): break

            element, mass, pseudo = line.split()
            self.settings_dict['pseudopotentials'].update({element : pseudo})
            i+=1

        #K_POINTS
        if 'gamma' in card_lines[k_points_index].split()[1].strip().lower():
            self.settings_dict['kpts'] = None
            self.settings_dict['koffset'] = None
        else:
            line = card_lines[k_points_index+1]
            self.settings_dict['kpts'] = list(map(int, line.split()[:3]))
            self.settings_dict['koffset'] = list(map(int, line.split()[3:]))


        self.settings_dict['additional_cards'] = []

        #HUBBARD
        if hubbard_index is not None:
            i = hubbard_index
            while i < len(card_lines):
                line = card_lines[i]

                if _end_of_card(card='HUBBARD', line=line): break

                self.settings_dict['additional_cards'] += [line]
                i+=1


@dataclass
class VaspSettings:
    '''
    Dataclass to store VASP input settings.
    '''
    vasp_pp_path: str
    vasp_pseudo_setups: Optional[dict]
    pymatgen_set: Optional[str]
    incar_path: Optional[str]
    kpoints_path: Optional[str]
    incar_path_screening: Optional[str]
    kpoints_path_screening: Optional[str]
    ediffg_screening: float = HYBRID_SCREENING_THRESHOLDS['vasp']
    vasp_xc_functional: str = "PBE"

    # non-initialized attributes (filled in __post_init__)
    incar_string: Optional[str] = None
    kpoints_string: Optional[str] = None
    incar_string_screening: Optional[str] = None
    kpoints_string_screening: Optional[str] = None
    settings_dict: Optional[dict] = None

    def __post_init__(self):

        #check parameters #############
        if self.pymatgen_set is None and self.incar_path is None and self.kpoints_path is None:
            raise RuntimeWarning('No pymatgen_set, incar_path or kpoints_path provided. '
            'Some default will be used, but they may not be suitable for your calculation.')

        if self.pymatgen_set is not None:
            self.pymatgen_set = self.pymatgen_set.lower()
            
            if self.pymatgen_set not in [
                'mprelaxset', 
                'mpmetalrelaxset', 
                'mpscanrelaxset', 
                'mphserelaxset', 
                'mitrelaxset']:
                raise ValueError(f'Invalid pymatgen_set: {self.pymatgen_set}. '
                                'Valid options are: MPRelaxSet, MPMetalRelaxSet, '
                                'MPScanRelaxSet, MPHSERelaxSet, MITRelaxSet (case insensitive).')


        #read incar and kpoints files if provided 
        if self.incar_path is not None:
            with open(self.incar_path, 'r',encoding=sys.getfilesystemencoding()) as f:
               self.incar_string = f.read()
        if self.kpoints_path is not None:
            with open(self.kpoints_path, 'r',encoding=sys.getfilesystemencoding()) as f:
                self.kpoints_string = f.read()
        if self.incar_path_screening is not None:
            with open(self.incar_path_screening, 'r',encoding=sys.getfilesystemencoding()) as f:
                self.incar_string_screening = f.read()
        if self.kpoints_path_screening is not None:
            with open(self.kpoints_path_screening, 'r',encoding=sys.getfilesystemencoding()) as f:
                self.kpoints_string_screening = f.read()

        #if any of them is present, create a settings_dict
        if self.incar_string or self.kpoints_string or self.incar_string_screening or self.kpoints_string_screening:
            self.settings_dict = {
                'incar': self.incar_string,
                'kpoints': self.kpoints_string,
                'incar_screening': self.incar_string_screening,
                'kpoints_screening': self.kpoints_string_screening,
            }
        

@dataclass
class MLSettings:
    '''
    Dataclass to store machine learning settings.
    '''
    force_conv_thr: float = 0.01


@dataclass
class DFTParams:
    '''
    Dataclass to store DFT program settings.
    '''

    program: str  # for convenience
    vasp: Optional[VaspSettings]
    espresso: Optional[EspressoSettings]
    ml: Optional[MLSettings]


    # non-initialized attributes (filled in __post_init__)
    settings_dict: Optional[dict] = None

    def __post_init__(self):
        self.program = self.program.lower()
        if self.program not in ['espresso', 'vasp', 'ml']:
            raise ValueError(f'DFT program must be one of ["espresso", "vasp", "ml"].')
        
        if self.program == 'espresso':
            if self.espresso is None:
                raise ValueError('espresso settings are missing.')
            else:
                self.settings_dict = self.espresso.settings_dict
        if self.program == 'vasp':
            if self.vasp is None:
                raise ValueError('vasp settings are missing.')
            else:
                self.settings_dict = self.vasp.settings_dict
        