#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#@author: Enrico Pedretti

"""
Module containing the Settings class,
used to read the settings file and store the parameters.

"""

from dataclasses import dataclass, field
from typing import Optional
from pathlib import Path
import tomllib
import json

from dacite import from_dict, Config

from xsorb.io.utils import ase_custom_read as read
from xsorb.dft_codes.definitions import SUPPORTED_PROGRAMS
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
    range_selection: dict

    atomic_species: Optional[list[str]]
    max_cn: Optional[float]
    max_cn_offset: float = 2
    include_surrounding_sites: bool = False
    surrounding_sites_deltaz: float = 1.5
    surrounding_exclude_main: bool = True
    cn_plain_fixed_radius: float = 1.5

    def __post_init__(self):
        if self.cn_method not in ['plain', 'minimumdistancenn', 'crystalnn']:
            raise ValueError('cn_method must be either plain, MinimumDistanceNN or CystalNN.')
        if 'mode' not in self.range_selection or 'value' not in self.range_selection:
            raise ValueError('range_selection must be in the format {mode = "max", value = 4}.')
        if self.range_selection['mode'] not in ['max', 'offset']:
            raise ValueError('range_selection mode must be either max or offset.')

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
    selected_atom_index: int
    x_rot_angles: list[float]
    y_rot_angles: list[float]
    z_rot_angles: list[float]

    individual_rotations: Optional[list[list[float]]]

    vertical_angles: str | list[float] = 'x'
    adsorption_distance_mode: str = 'value'
    target_distance: float = 2.0
    min_distance: float = 1.5
    radius_scale_factor: float = 1.0

    def __post_init__(self):
        if 'mode' not in self.molecule_axis or 'values' not in self.molecule_axis:
            raise ValueError('molecule_axis must be in the format \
                             {mode = "atom_indices", values = [1, 2]} or\
                              {mode = "vector", values = [1,0,0]}.')
        if self.molecule_axis['mode'] not in ['atom_indices', 'vector']:
            raise ValueError('molecule_axis mode must be either atom_indices or vector.')
        if self.molecule_axis['mode'] == 'atom_indices':
            if len(self.molecule_axis['values']) != 2:
                raise ValueError('molecule_axis values must be a list of two atom indices \
                                 when mode is atom_indices.')
        elif self.molecule_axis['mode'] == 'vector':
            if len(self.molecule_axis['values']) != 3:
                raise ValueError('molecule_axis values must be a list of three floats \
                                 when mode is vector.')

        if self.adsorption_distance_mode is not None:
            if self.adsorption_distance_mode not in ['value', 'covalent_radius', 'vdw_radius']:
                raise ValueError('adsorption_distance_mode must be either value, \
                                 covalent_radius or vdw_radius.')

        if self.target_distance < self.min_distance:
            raise ValueError('target_distance must be greater than min_distance.')

        if isinstance(self.vertical_angles, list):
            if len(self.vertical_angles) == 0:
                raise ValueError('vertical_angles given as a list has to contain at least 1 angle.')

        elif isinstance(self.vertical_angles, str):
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
    fix_slab_xyz: list[int] = field(default_factory=lambda: [0, 0, 0])
    fix_mol_xyz: list[int] = field(default_factory=lambda: [0, 0, 1])
    fix_bondlengths_preopt: bool = False
    fix_slab_preopt: bool = False

    def __post_init__(self):
        if self.fixed_layers_slab is not None and self.fixed_indices_slab is not None:
            raise ValueError('You can use either fixed_layers_slab or fixed_indices_slab, \
                             not both at the same time.')
        if self.fix_bondlengths_preopt and self.fix_slab_preopt:
            raise ValueError('fix_bondlengths_preopt and fix_slab_preopt cannot be used together.')

@dataclass
class MiscParams:
    '''
    Dataclass to store the parameters in the
    misc card of the settings file.
    '''
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
    misc: MiscParams = field(default_factory=MiscParams())



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
            with open("settings.json", "rb") as f:
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
                               data=self.settings_dict["Input"],
                               config=Config(type_hooks={str: str.lower}, strict=True))


        self.structure = from_dict(data_class=StructureParams,
                                   data=self.settings_dict["Structure"],
                                   config=Config(type_hooks={str: str.lower}, strict=True))


        #initialize the dft settings
        self.program : str = self.settings_dict["Calculation_parameters"]["dft_program"].upper()

        if self.program not in SUPPORTED_PROGRAMS:
            raise ValueError('dft_program must be either "VASP" or "ESPRESSO".')
        if self.program not in self.settings_dict["Calculation_parameters"]:
            raise ValueError(f"Settings for {self.program} are missing.")

        self.dftprogram_settings_dict = get_dftprogram_settings(
            self.program,
            self.settings_dict["Calculation_parameters"][self.program]
            )

        #read energies of slab and mol if requested
        if read_energies and self.input.E_slab_mol is None:
            self.read_E_slab_mol(verbose)


    def read_E_slab_mol(self, verbose : bool = True): # pylint: disable=invalid-name
        '''
        Attempt to read the energies of the slab and molecule from the slab and molecule files,
        and store them in E_slab_mol of the input dataclass.
        If the energies cannot be read, E_slab_mol is set to [0.0, 0.0].
        '''
        try:
            slab_en = read(self.input.slab_filename).get_potential_energy()
            mol_en = read(self.input.molecule_filename).get_potential_energy()
            self.input.E_slab_mol = [slab_en, mol_en]
        except Exception as e: # pylint: disable=broad-except
            #this is a general exception to catch any error that might occur.
            #It can be file not found, or energy not present in the file, etc.
            self.input.E_slab_mol = [0.0, 0.0]
            if verbose:
                print('It was not possible to obtain slab and molecule energy in any way.',
                       f'Error message from ase: {e}.',
                        'Total energies will be shown instead of adsorption energies.')
