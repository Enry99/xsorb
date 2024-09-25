#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#@author: Enrico Pedretti

"""
Class to store the parameters for generating the geometries and launching the calculations.

"""

from dataclasses import dataclass, field
from typing import Optional
from dacite import from_dict, Config
import tomllib
from xsorb.dft_codes.definitions import SUPPORTED_PROGRAMS
from xsorb.dft_codes.input_settings import get_dftprogram_settings


@dataclass
class Input:
    slab_filename: str
    molecule_filename: str
    jobscript_path: str
    submit_command: str
    E_slab_mol : Optional[list[float]]
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
class High_symmetry_params:
    symm_reduce: float = 0.01
    near_reduce: float = 0.01

@dataclass
class Coord_number_params:
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
class Adsorption_sites:
    mode : str    
    selected_sites: Optional[list[int]]
    high_symmetry_params: Optional[High_symmetry_params]
    coord_number_params: Optional[Coord_number_params]
    surface_height: float = 0.9

    def __post_init__(self):
        if self.mode not in ['high_symmetry', 'coord_number']:
            raise ValueError('mode must be either high_symmetry or coord_number.')
        if self.mode == 'high_symmetry':
            if self.high_symmetry_params is None:
                raise ValueError('high_symmetry_params must be provided when mode is high_symmetry.')
        elif self.mode == 'coord_number':
            if self.coord_number_params is None:
                raise ValueError('coord_number_params must be provided when mode is coord_number.')

@dataclass
class Molecule:
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
    radius_scale_factor: float = 1.1

    def __post_init__(self):
        if 'mode' not in self.molecule_axis or 'values' not in self.molecule_axis:
            raise ValueError('molecule_axis must be in the format {mode = "atom_indices", values = [1, 2]} or {mode = "vector", values = [1,0,0]}.')
        if self.molecule_axis['mode'] not in ['atom_indices', 'vector']:
            raise ValueError('molecule_axis mode must be either atom_indices or vector.')
        if self.molecule_axis['mode'] == 'atom_indices':
            if len(self.molecule_axis['values']) != 2:
                raise ValueError('molecule_axis values must be a list of two atom indices when mode is atom_indices.')
        elif self.molecule_axis['mode'] == 'vector':
            if len(self.molecule_axis['values']) != 3:
                raise ValueError('molecule_axis values must be a list of three floats when mode is vector.')
            
        if self.adsorption_distance_mode is not None:
            if self.adsorption_distance_mode not in ['value', 'covalent_radius', 'vdw_radius']:
                raise ValueError('adsorption_distance_mode must be either value, covalent_radius or vdw_radius.')

        if type(self.vertical_angles) is list:
            if len(self.vertical_angles) == 0:
                raise ValueError('vertical_angles given as a list has to contain at least one angle.')

        elif type(self.vertical_angles) is str:
            if self.vertical_angles not in ['x', 'z', 'none']:
                raise ValueError('vertical_angles, when not given as a list, must be either "x", "z", "none".')
            
            if self.vertical_angles == 'x': self.vertical_angles = self.x_rot_angles
            elif self.vertical_angles == 'z': self.vertical_angles = self.z_rot_angles
            elif self.vertical_angles == 'none': self.vertical_angles = None
                

            
@dataclass
class Constraints:
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
            raise ValueError('You can use either fixed_layers_slab or fixed_indices_slab, not both at the same time.')
        if self.fix_bondlengths_preopt and self.fix_slab_preopt:
            raise ValueError('fix_bondlengths_preopt and fix_slab_preopt cannot be used together.')
        
@dataclass
class Misc:
    mol_before_slab: bool = False
    sort_atoms_by_z: bool = True
    translate_slab: bool = True

@dataclass
class Structure:
    adsorption_sites: Adsorption_sites
    molecule: Molecule
    constraints: Constraints = field(default_factory=lambda: Constraints(None, None, None))
    misc: Misc = field(default_factory=lambda: Misc())



class Settings:
    
    def __init__(self, 
                 settings_filename: str = "settings.toml", 
                 READ_ENERGIES: bool = True, 
                 VERBOSE: bool = True):
        
        # read the settings file
        with open(settings_filename, "rb") as f:
            self.settings_dict = tomllib.load(f)

        #check for existence of the main cards
        cards = ['Input','Structure', 'Calculation_parameters']
        for card in cards:
            if card not in self.settings_dict:
                raise RuntimeError(f"{card} card not found in {settings_filename}.")    


        #intialize the dataclasses
        self.input = from_dict(data_class=Input, 
                               data=self.settings_dict["Input"], 
                               config=Config(type_hooks={str: str.lower}, strict=True))
    
        
        self.structure = from_dict(data_class=Structure, 
                                   data=self.settings_dict["Structure"], 
                                   config=Config(type_hooks={str: str.lower}, strict=True))
    

        #initialize the dft settings
        self.program : str = self.settings_dict["Calculation_parameters"]["dft_program"].upper()
        
        if self.program not in SUPPORTED_PROGRAMS:
            raise ValueError('dft_program must be either "VASP" or "ESPRESSO".')
        if self.program not in self.settings_dict["Calculation_parameters"]:
            raise ValueError(f"Settings for {self.program} are missing.")

        self.dftprogram_settings_dict = get_dftprogram_settings(self.program, self.settings_dict["Calculation_parameters"][self.program])
        
        #read energies of slab and mol if requested
        if READ_ENERGIES and self.input.E_slab_mol is None:
            self.read_E_slab_mol(VERBOSE)

                    
    def read_E_slab_mol(self, VERBOSE : bool = True):
        try:
            from ase.io import read
            import ase_custom
            slab_en = read(self.input.slab_filename).get_potential_energy()
            mol_en = read(self.input.molecule_filename).get_potential_energy()
            self.input.E_slab_mol = [slab_en, mol_en]
        except:
            self.input.E_slab_mol = [0.0, 0.0]
            if VERBOSE:
                print('It was not possible to obtain slab and molecule energy in any way.' \
                        'Total energies will be shown instead of adsorption energies.')
