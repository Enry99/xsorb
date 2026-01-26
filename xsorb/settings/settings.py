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
from glob import glob
from dataclasses import dataclass

try:
    import tomllib
except ModuleNotFoundError:
    import pip._vendor.tomli as tomllib

from dacite import from_dict, Config

from xsorb.settings.input_settings import InputParams, StructureParams, DatabaseParams
from xsorb.dft_codes.input_settings import DFTParams
from xsorb.dft_codes.definitions import OUT_FILE_PATHS
from xsorb.ase_custom.io import ase_custom_read as read


@dataclass
class Energies:
    '''
    Dataclass to store the "final" energies of the slab and molecule
    after reading settings and calculation files.
    Values are set to 0 if no actual values are provided.
    '''
    E_slab : float = 0.0 # pylint: disable=invalid-name
    E_mol : float|list[float] = 0.0 # pylint: disable=invalid-name
    E_slab_ml : float = 0.0 # pylint: disable=invalid-name
    E_mol_ml : float|list[float] = 0.0 # pylint: disable=invalid-name



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

    # Attibutes:
    # input: InputParams
    # structure: StructureParams
    # dft: DFTParams
    # database: DatabaseParams
    # energies: Energies

    def __init__(self,
                 read_energies_dft: bool = False,
                 read_energies_ml: bool = False,
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

        self.database = from_dict(data_class=DatabaseParams,
                                    data=settings_dict.get("Database", {}),)

        self.energies = Energies()

        ################################
        #read energies if requested
        if read_energies_dft:
            self.read_E_slab_mol(verbose)
        if read_energies_ml:
            self.read_E_slab_mol_ml(verbose)


    def read_E_slab_mol(self, verbose : bool = True): # pylint: disable=invalid-name
        '''
        Attempt to read the energies of the slab and molecule from the slab and molecule files,
        and store them in the input dataclass.
        If either file is not found, the corresponding energy will be set to 0.0.
        '''

        if self.input.E_slab is None:
            try:
                self.energies.E_slab = \
                    read(OUT_FILE_PATHS['slab'][self.dft.program]).get_potential_energy()
            except Exception: # pylint: disable=broad-except
                try:
                    self.energies.E_slab = read(self.input.slab_filename).get_potential_energy()
                except Exception as e: # pylint: disable=broad-except
                    if verbose:
                        logging.error(f"Error reading slab energy: {e}. Setting to 0") # pylint: disable=logging-fstring-interpolation
        if self.input.E_mol is None:
            try:
                n_conformers = len(glob(OUT_FILE_PATHS['mol'][self.dft.program].format('*')))
                energies = [read(
                    OUT_FILE_PATHS['mol'][self.dft.program].format(i)).get_potential_energy()
                    for i in range(n_conformers)]
                self.energies.E_mol = energies[0] if len(energies) == 1 else energies
            except Exception: # pylint: disable=broad-except
                try:
                    traj = read(self.input.molecule_filename, index=':')
                    energies = [at.get_potential_energy() for at in traj]
                    self.energies.E_mol = energies[0] if len(energies) == 1 else energies
                except Exception as e: # pylint: disable=broad-except
                    if verbose:
                        logging.error(f"Error reading molecule energy: {e}. Setting to 0") # pylint: disable=logging-fstring-interpolation


    def read_E_slab_mol_ml(self, verbose : bool = True): # pylint: disable=invalid-name
        '''
        Attempt to read the energies of the slab and molecule from the slab and molecule files,
        and store them in the input dataclass.
        If either file is not found, the corresponding energy will be set to 0.0.
        '''

        try:
            self.energies.E_slab_ml = read(OUT_FILE_PATHS['slab']['ml']).get_potential_energy()
        except Exception as e: # pylint: disable=broad-except
            if verbose:
                logging.error(f"Error reading ML slab energy: {e}. Setting to 0") # pylint: disable=logging-fstring-interpolation
        try:
            n_conformers = len(glob(OUT_FILE_PATHS['mol']['ml'].format('*')))
            #read preserving order
            energies = [read(OUT_FILE_PATHS['mol']['ml'].format(i)).get_potential_energy()
                            for i in range(n_conformers)]
            self.energies.E_mol_ml = energies[0] if len(energies) == 1 else energies
        except Exception as e: # pylint: disable=broad-except
            if verbose:
                logging.error(f"Error reading ML molecule energy: {e}. Setting to 0") # pylint: disable=logging-fstring-interpolation
