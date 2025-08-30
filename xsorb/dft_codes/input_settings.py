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
from copy import deepcopy

from ase.io.espresso import read_fortran_namelist

from xsorb.dft_codes.definitions import HYBRID_SCREENING_THRESHOLDS

# dataclasses to store DFT program settings
@dataclass
class EspressoParams:
    '''
    Dataclass to store Espresso input settings.
    '''
    pwi_path: str
    pwi_path_screening : Optional[str]
    etot_conv_thr_screening: float = HYBRID_SCREENING_THRESHOLDS['espresso'][0]
    forc_conv_thr_screening: float = HYBRID_SCREENING_THRESHOLDS['espresso'][1]

    # non-initialized attributes (filled in __post_init__)
    settings_dict = None
    settings_dict_screening = None

    def build_settings_dict(self, filepath : str):
        '''
        Build the settings dictionary from the input files.
        '''
        #NOTE 1: The blocks CELL_PARAMETERS ATOMIC_POSITIONS ATOMIC_SPECIES must NOT be included
        # in input file, as they are read from the input structures
        #NOTE 2: This code does not yet support the following Espresso blocks:
        #OCCUPATIONS, CONSTRAINTS, ATOMIC_VELOCITIES, ATOMIC_FORCES, ADDITIONAL_K_POINTS, SOLVENTS

        # parse namelist section and extract remaining lines
        with open(filepath, 'r') as file:
            _settings_dict, card_lines = read_fortran_namelist(file)
            _settings_dict = dict(_settings_dict)


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
        _settings_dict['pseudopotentials'] = {}
        i = atomic_species_index+1
        while i < len(card_lines):
            line = card_lines[i]

            if _end_of_card(card='ATOMIC_SPECIES', line=line): break

            element, mass, pseudo = line.split()
            _settings_dict['pseudopotentials'].update({element : pseudo})
            i+=1

        #K_POINTS
        if 'gamma' in card_lines[k_points_index].split()[1].strip().lower():
            _settings_dict['kpts'] = None
            _settings_dict['koffset'] = None
        else:
            line = card_lines[k_points_index+1]
            _settings_dict['kpts'] = list(map(int, line.split()[:3]))
            _settings_dict['koffset'] = list(map(int, line.split()[3:]))


        _settings_dict['additional_cards'] = []

        #HUBBARD
        if hubbard_index is not None:
            i = hubbard_index
            while i < len(card_lines):
                line = card_lines[i]

                if _end_of_card(card='HUBBARD', line=line): break

                _settings_dict['additional_cards'] += [line]
                i+=1

        return _settings_dict


    def __post_init__(self):

        #read pwi file and build settings_dict
        self.settings_dict = self.build_settings_dict(self.pwi_path)
        self.settings_dict['control'].update({'calculation': 'relax' })
        self.settings_dict['control'].update({'restart_mode': 'from_scratch'})
        self.settings_dict['control'].update({'outdir': 'OUT'})
        if 'ions' not in self.settings_dict: self.settings_dict['ions'] = {}

        #if screening pwi file is provided, build screening settings_dict
        if self.pwi_path_screening is not None:
            self.settings_dict_screening = self.build_settings_dict(self.pwi_path_screening)
            self.settings_dict_screening['control'].update({'calculation': 'relax' })
            self.settings_dict_screening['control'].update({'restart_mode': 'from_scratch'})
            self.settings_dict_screening['control'].update({'outdir': 'OUT'})
            if 'ions' not in self.settings_dict_screening: self.settings_dict['ions'] = {}
        else:
            self.settings_dict_screening = deepcopy(self.settings_dict)
        #update screening settings_dict with etot_conv_thr and forc_conv_thr
        self.settings_dict_screening['control'].update({'etot_conv_thr': self.etot_conv_thr_screening})
        self.settings_dict_screening['control'].update({'forc_conv_thr': self.forc_conv_thr_screening})


@dataclass
class VaspParams:
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
    settings_dict = None
    settings_dict_screening = None

    def build_settings_dict(self, incar_path: str | None, kpoints_path: str | None):
        '''
        Build the settings dictionary from the input files.
        '''
        #parse INCAR and KPOINTS files
        _settings_dict = {}
        if incar_path is not None:
            with open(incar_path, 'r',encoding=sys.getfilesystemencoding()) as f:
               incar_string = f.read()
               _settings_dict['incar_string'] = incar_string

        if kpoints_path is not None:
            with open(kpoints_path, 'r',encoding=sys.getfilesystemencoding()) as f:
                kpoints_string = f.read()
                _settings_dict['kpoints_string'] = kpoints_string

        #add vasp_pp_path, vasp_pseudo_setups and pymatgen_set and vasp_xc_functional
        # to the settings dictionary
        _settings_dict['vasp_pp_path'] = self.vasp_pp_path
        _settings_dict['vasp_pseudo_setups'] = self.vasp_pseudo_setups
        _settings_dict['pymatgen_set'] = self.pymatgen_set
        _settings_dict['vasp_xc_functional'] = self.vasp_xc_functional

        return _settings_dict


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
        if self.incar_path is not None or self.kpoints_path is not None:
            self.settings_dict = self.build_settings_dict(self.incar_path, self.kpoints_path)

            #fix if user forgot to put the correct IBRION for relax
            if 'incar_string' in self.settings_dict and self.pymatgen_set is None:
                s = self.settings_dict['incar_string'].split('\n')

                missing_ibrion = True
                for i, line in enumerate(s):
                    if 'IBRION' in line: missing_ibrion = False
                if missing_ibrion: s.append('IBRION = 2')

                self.settings_dict['incar_string'] = '\n'.join(s)


        if self.incar_path_screening is not None or self.kpoints_path_screening is not None:
            self.settings_dict_screening = self.build_settings_dict(
                self.incar_path_screening, self.kpoints_path_screening)

            #fix if user forgot to put the correct IBRION for relax, and add EDIFFG for screening.
            #EDIFFG is put here so that in any case this will be the final value, even if the user had
            #specified a different value explicitly in the INCAR instead of using one of the RelaxSets
            if 'incar_string' in self.settings_dict_screening and self.pymatgen_set is None:
                s = self.settings_dict_screening['incar_string'].split('\n')
            else: s = []

            missing_ibrion = True
            missing_ediffg = True
            for i, line in enumerate(s):
                if 'EDIFFG' in line:
                    s[i] = f"EDIFFG = {self.ediffg_screening}"
                    missing_ediffg = False
                if 'IBRION' in line:
                    missing_ibrion = False
            if missing_ediffg: s.append(f"EDIFFG = {self.ediffg_screening}")
            if missing_ibrion: s.append('IBRION = 2')

            self.settings_dict_screening['incar_string'] = '\n'.join(s)


@dataclass
class MLParams:
    '''
    Dataclass to store machine learning settings.
    '''
    force_conv_thr: float = 0.01 # force threshold (in eV/A). NOT IMPLEMENTED YET!


@dataclass
class DFTParams:
    '''
    Dataclass to store DFT program settings.
    '''

    program: str  # for convenience
    vasp: Optional[VaspParams]
    espresso: Optional[EspressoParams]
    ml: Optional[MLParams]


    def __post_init__(self):
        self.program = self.program.lower()
        if self.program not in ['espresso', 'vasp', 'ml']:
            raise ValueError(f'DFT program must be one of ["espresso", "vasp", "ml"].')

        if self.program == 'espresso':
            if self.espresso is None:
                raise ValueError('espresso settings are missing.')
        if self.program == 'vasp':
            if self.vasp is None:
                raise ValueError('vasp settings are missing.')


    def get_settings_dict(self, calc_type: str = 'relax', force_gamma : bool = False) -> dict:
        '''
        Get the settings dictionary for the specified calculation type.
        '''
        if self.program == 'espresso':
            if calc_type == 'screening' and self.espresso.settings_dict_screening is not None:
                settings_dict = deepcopy(self.espresso.settings_dict_screening)
            else:
                settings_dict = deepcopy(self.espresso.settings_dict)

            if force_gamma:
                settings_dict['kpts'] = None
                settings_dict['koffset'] = None

            return settings_dict

        elif self.program == 'vasp':
            if calc_type == 'screening' and self.vasp.settings_dict_screening is not None:
                settings_dict = deepcopy(self.vasp.settings_dict_screening)
            else:
                settings_dict = deepcopy(self.vasp.settings_dict)

            if force_gamma:
                settings_dict['kpoints_string'] = 'Gamma-point only\n0\nMonkhorst Pack\n1 1 1\n0 0 0'

            return settings_dict

        elif self.program == 'ml':
            return {}
            #return {'force_conv_thr': self.ml.force_conv_thr}

        else:
            raise ValueError(f'Program {self.program} not recognized.')
