#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#@author: Enrico Pedretti

"""
Module containing the Settings class,
used to read the settings file and store the parameters.

"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional, Union
from pathlib import Path
import json
import sys

try:
    import tomllib
except ModuleNotFoundError:
    import pip._vendor.tomli as tomllib

from dacite import from_dict, Config

from xsorb.ase_custom.io import ase_custom_read as read
from xsorb.dft_codes.definitions import SUPPORTED_PROGRAMS, OUT_FILE_PATHS
from xsorb.dft_codes.input_settings import get_dftprogram_settings


@dataclass
class InputParams:
    '''
    Dataclass to store the parameters in the
    INPUT card of the settings file.
    '''

    slab_filename: str
    molecule_filename: str
    jobscript_path: str
    submit_command: str
    E_slab_mol : Optional[list[float]] # pylint: disable=invalid-name
    jobscript_ml_path: Optional[str]
    submit_command_ml: Optional[str]
    jobname_prefix: str = ''

    def __post_init__(self):
        if self.E_slab_mol is not None:
            if len(self.E_slab_mol) != 2:
                raise ValueError("E_slab_mol must be a list of two floats.")

        if self.jobscript_ml_path is not None and self.submit_command_ml is None:
            raise ValueError("jobscript_ml_path is provided but submit_command_ml is not.")

@dataclass
class HighSymmetryParams:
    '''
    Dataclass to store the parameters in the
    high_symmetry_params card of the settings file.
    '''
    symm_reduce: float = 0.01

@dataclass
class CoordNumberParams:
    '''
    Dataclass to store the parameters in the
    coord_number_params card of the settings file.
    '''
    cn_method: str

    atomic_species: Optional[list[str]]
    max_cn: Optional[float]
    max_cn_offset: Optional[float]
    include_surrounding_sites: bool = False
    surrounding_sites_deltaz: float = 1.5
    cn_plain_fixed_radius: float = 1.5

    def __post_init__(self):
        if self.cn_method not in ['plain', 'minimumdistancenn', 'crystalnn']:
            raise ValueError('cn_method must be either plain, MinimumDistanceNN or CystalNN.')
        if self.max_cn is not None and self.max_cn_offset is not None:
            raise ValueError('You can use either max_cn or max_cn_offset, not both at the same time.')
        if self.max_cn is None and self.max_cn_offset is None:
            self.max_cn_offset = 2

@dataclass
class AdsorptionSitesParams:
    '''
    Dataclass to store the parameters in the
    adsorption_sites card of the settings file.
    '''
    mode : str
    selected_sites: Optional[list[int]]
    high_symmetry_params: Optional[HighSymmetryParams]
    coord_number_params: Optional[CoordNumberParams]
    surface_thickness: float = 0.9

    def __post_init__(self):
        if self.mode not in ['high_symmetry', 'coord_number']:
            raise ValueError('mode must be either high_symmetry or coord_number.')
        if self.mode == 'high_symmetry':
            if self.high_symmetry_params is None:
                raise ValueError('high_symmetry_params must be provided when mode is high_symmetry')
        elif self.mode == 'coord_number':
            if self.coord_number_params is None:
                raise ValueError('coord_number_params must be provided when mode is coord_number.')

@dataclass
class MoleculeParams:
    '''
    Dataclass to store the parameters in the
    molecule card of the settings file.
    '''
    molecule_axis: dict
    selected_atom_indexes: list[int]
    x_rot_angles: list[float]
    y_rot_angles: list[float]
    z_rot_angles: list[float]

    individual_rotations: Optional[list[list[float]]]

    vertical_angles: Union[str,list[float]] = 'x'
    adsorption_distance_mode: str = 'value'
    target_distance: float = 2.0
    min_distance: float = 1.5
    radius_scale_factor: float = 1.1

    def __post_init__(self):
        if 'mode' not in self.molecule_axis or 'values' not in self.molecule_axis:
            raise ValueError('molecule_axis must be in the format \
                             {mode = "atom_indices", values = [1, 2]} or\
                              {mode = "vector", values = [1,0,0]}.')
        if self.molecule_axis['mode'] not in ['atom_indices', 'vector']:
            raise ValueError('molecule_axis mode must be either atom_indices or vector.')
        if self.molecule_axis['mode'] == 'atom_indices' and len(self.molecule_axis['values']) != 2:
            raise ValueError('molecule_axis values must be a list of two atom indices \
                                when mode is atom_indices.')
        if self.molecule_axis['mode'] == 'vector' and len(self.molecule_axis['values']) != 3:
            raise ValueError('molecule_axis values must be a list of three floats \
                                when mode is vector.')

        if self.adsorption_distance_mode is not None and \
            self.adsorption_distance_mode not in ['value', 'covalent_radius', 'vdw_radius']:
            raise ValueError('adsorption_distance_mode must be either value, \
                                covalent_radius or vdw_radius.')

        if self.target_distance < self.min_distance:
            raise ValueError('target_distance must be greater than min_distance.')

        if isinstance(self.vertical_angles, list) and len(self.vertical_angles) == 0:
            raise ValueError('vertical_angles given as a list has to contain at least 1 angle.')

        if isinstance(self.vertical_angles, str):
            if self.vertical_angles not in ['x', 'z', 'none']:
                raise ValueError('vertical_angles, when not given as a list, must be either \
                                 "x", "z", "none".')

            if self.vertical_angles == 'x':
                self.vertical_angles = self.x_rot_angles
            elif self.vertical_angles == 'z':
                self.vertical_angles = self.z_rot_angles
            elif self.vertical_angles == 'none':
                self.vertical_angles = None



@dataclass
class ConstraintsParams:
    '''
    Dataclass to store the parameters in the
    constraints card of the settings file.
    '''
    fixed_layers_slab: Optional[list[int]]
    fixed_indices_slab: Optional[list[int]]
    fixed_indices_mol: Optional[list[int]]
    layers_height: float = 0.5
    fix_slab_xyz: list[bool] = field(default_factory=lambda: [True,True,True])
    fix_mol_xyz: list[bool] = field(default_factory=lambda: [True,True,False])
    fix_slab_ml_opt: bool = False

    def __post_init__(self):
        if self.fixed_layers_slab is not None and self.fixed_indices_slab is not None:
            raise ValueError('You can use either fixed_layers_slab or fixed_indices_slab, \
                             not both at the same time.')

@dataclass
class MiscParams:
    '''
    Dataclass to store the parameters in the
    misc card of the settings file.
    '''
    inside_only: bool = False
    mol_before_slab: bool = False
    sort_atoms_by_z: bool = True
    translate_slab: bool = True

@dataclass
class StructureParams:
    '''
    Dataclass to store the parameters in the
    structure card of the settings file.
    '''
    adsorption_sites: AdsorptionSitesParams
    molecule: MoleculeParams
    constraints: ConstraintsParams = field(
        default_factory=lambda: ConstraintsParams(None, None, None))
    misc: MiscParams = field(default_factory=lambda: MiscParams(False,False,True,True))



class Settings:
    '''
    Class to read the settings file and store all the input parameters.
    Checks are performed to ensure that the settings file is correctly formatted.

    The file must be present in the working directory,
    either in .toml or .json format: settings.toml or settings.json.

    Initialization parameters:
    - read_energies : if True, it will attempt to read energies of the slab and molecule
    - verbose : if True, messages will be printed to standard output
    '''

    def __init__(self,
                 read_energies: bool = False,
                 verbose: bool = True):

        #Read the settings file
        if Path('settings.toml').is_file():
            with open("settings.toml", "rb") as f:
                self.settings_dict = tomllib.load(f)
        elif Path('settings.json').is_file():
            with open("settings.json", "r", encoding=sys.getfilesystemencoding()) as f:
                self.settings_dict = json.load(f)
        else:
            raise FileNotFoundError("Settings file (settings.toml or settings.json)"\
                                    " not found in working directory. Terminating.")

        #check for existence of the main cards
        cards = ['Input','Structure', 'Calculation_parameters']
        for card in cards:
            if card not in self.settings_dict:
                raise RuntimeError(f"{card} card not found in settings file.")


        #intialize the dataclasses
        self.input = from_dict(data_class=InputParams,
                               data=self.settings_dict["Input"])


        self.structure = from_dict(data_class=StructureParams,
                                   data=self.settings_dict["Structure"],
                                   config=Config(type_hooks={str: str.lower}, strict=True))


        #initialize the dft settings
        self.program : str = self.settings_dict["Calculation_parameters"]["dft_program"].lower()

        if self.program not in SUPPORTED_PROGRAMS:
            raise ValueError(f'DFT program must be one of {SUPPORTED_PROGRAMS}.')
        if self.program not in self.settings_dict["Calculation_parameters"]:
            raise ValueError(f"Settings for {self.program} are missing.")

        self.dftprogram_settings_dict = get_dftprogram_settings(
            self.program,
            self.settings_dict["Calculation_parameters"][self.program]
            )

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
                self.input.E_slab_mol[0] = read(self.input.slab_filename).get_potential_energy()
            except Exception as e: # pylint: disable=broad-except
                try:
                    self.input.E_slab_mol[0] = \
                        read(OUT_FILE_PATHS['slab'][self.program]).get_potential_energy()
                except Exception as e: # pylint: disable=broad-except
                    if verbose:
                        print(f"Error reading slab energy: {e}. Setting to 0")
        if int(self.input.E_slab_mol[1]) == 0:
            try:
                self.input.E_slab_mol[1] = read(self.input.molecule_filename).get_potential_energy()
            except Exception as e: # pylint: disable=broad-except
                try:
                    self.input.E_slab_mol[1] = \
                        read(OUT_FILE_PATHS['mol'][self.program]).get_potential_energy()
                except Exception as e: # pylint: disable=broad-except
                    if verbose:
                        print(f"Error reading molecule energy: {e}. Setting to 0")

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
                print(f"Error reading ML slab energy: {e}. Setting to 0")
        try:
            emol_ml = \
                read(OUT_FILE_PATHS['mol']['ml']).get_potential_energy()
        except Exception as e: # pylint: disable=broad-except
            emol_ml = 0.0
            if verbose:
                print(f"Error reading molecule energy: {e}. Setting to 0")

        return eslab_ml + emol_ml
