#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#@author: Enrico Pedretti

"""
Module containing the Settings class,
used to read the settings file and store the parameters.

"""

from __future__ import annotations
from pathlib import Path
import json
import sys
import logging

try:
    import tomllib
except ModuleNotFoundError:
    import pip._vendor.tomli as tomllib

from dacite import from_dict, Config

from xsorb.settings.input_settings import InputParams, StructureParams
from xsorb.dft_codes.input_settings import DFTParams
from xsorb.dft_codes.definitions import OUT_FILE_PATHS
from xsorb.ase_custom.io import ase_custom_read as read


class Settings:
    '''
    Class to read the settings file and store all the input parameters.
    Checks are performed to ensure that the settings file is correctly formatted.

    The file must be present in the working directory,
    either in .toml or .json format: settings.toml or settings.json.

    Initialization parameters:
    - read_energies : if True, it will attempt to read energies of the slab and molecule
    - verbose : if True, messages will be printed to standard output

    Attributes:
    - input : InputParams dataclass containing the input parameters
    - structure : StructureParams dataclass containing the structure parameters
    - dft : DFTParams dataclass containing the DFT program settings
    '''

    input: InputParams
    structure: StructureParams
    dft: DFTParams

    def __init__(self,
                 read_energies: bool = False,
                 verbose: bool = True):

        #Read the settings file
        if Path('settings.toml').is_file():
            with open("settings.toml", "rb") as f:
                settings_dict = tomllib.load(f)
        elif Path('settings.json').is_file():
            with open("settings.json", "r", encoding=sys.getfilesystemencoding()) as f:
                settings_dict = json.load(f)
        else:
            raise FileNotFoundError("Settings file (settings.toml or settings.json)"\
                                    " not found in working directory. Quitting.")

        ################################
        #check for existence of the main cards
        cards = ['Input','Structure', 'Calculation_parameters']
        for card in cards:
            if card not in settings_dict:
                raise RuntimeError(f"{card} card not found in settings file.")

        ################################
        #intialize the dataclasses
        self.input = from_dict(data_class=InputParams,
                               data=settings_dict["Input"])

        self.structure = from_dict(data_class=StructureParams,
                                   data=settings_dict["Structure"],
                                   config=Config(type_hooks={str: str.lower}, strict=True))

        self.dft = from_dict(data_class=DFTParams,
                             data=settings_dict["Calculation_parameters"])

        ################################
        #at this point, self.input.E_slab_mol can be None, if not specified in the settings file,
        #or a list of two floats, if specified. One can be 0, e.g. [13.6, 0.0] if only the
        #slab energy or molecule energy is known. We need to fill in the missing energies if
        #available if read_energies is True.
        if read_energies:
            self.total_e_slab_mol = self.read_E_slab_mol(verbose)
            self.total_e_slab_mol_ml = self.read_E_slab_mol_ml(verbose)
        else:
            self.total_e_slab_mol = None
            self.total_e_slab_mol_ml = None


    def read_E_slab_mol(self, verbose : bool = True): # pylint: disable=invalid-name
        '''
        Attempt to read the energies of the slab and molecule from the slab and molecule files,
        and store them in E_slab_mol of the input dataclass.
        If either file is not found, the corresponding energy will be set to 0.0.
        '''

        if self.input.E_slab_mol is None:
            self.input.E_slab_mol = [0,0]

        if int(self.input.E_slab_mol[0]) == 0:
            try:
                self.input.E_slab_mol[0] = \
                    read(OUT_FILE_PATHS['slab'][self.dft.program]).get_potential_energy()
            except Exception: # pylint: disable=broad-except
                try:
                    self.input.E_slab_mol[0] = read(self.input.slab_filename).get_potential_energy()
                except Exception as e: # pylint: disable=broad-except
                    if verbose:
                        logging.error(f"Error reading slab energy: {e}. Setting to 0")
        if int(self.input.E_slab_mol[1]) == 0:
            try:
                self.input.E_slab_mol[1] = \
                    read(OUT_FILE_PATHS['mol'][self.dft.program]).get_potential_energy()

            except Exception: # pylint: disable=broad-except
                try:
                    self.input.E_slab_mol[1] = read(self.input.molecule_filename).get_potential_energy()
                except Exception as e: # pylint: disable=broad-except
                    if verbose:
                        logging.error(f"Error reading molecule energy: {e}. Setting to 0")

        return sum(self.input.E_slab_mol)


    def read_E_slab_mol_ml(self, verbose : bool = True): # pylint: disable=invalid-name
        '''
        Attempt to read the energies of the slab and molecule from the slab and molecule files,
        and store them in E_slab_mol of the input dataclass.
        If either file is not found, the corresponding energy will be set to 0.0.
        '''

        try:
            eslab_ml = \
                read(OUT_FILE_PATHS['slab']['ml']).get_potential_energy()
        except Exception as e: # pylint: disable=broad-except
            eslab_ml = 0.0
            if verbose:
                logging.error(f"Error reading ML slab energy: {e}. Setting to 0")
        try:
            emol_ml = \
                read(OUT_FILE_PATHS['mol']['ml']).get_potential_energy()
        except Exception as e: # pylint: disable=broad-except
            emol_ml = 0.0
            if verbose:
                logging.error(f"Error reading molecule energy: {e}. Setting to 0")

        return eslab_ml + emol_ml
