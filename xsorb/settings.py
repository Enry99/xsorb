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
#import xsorb
from dftcode_specific import SUPPORTED_PROGRAMS, HYBRID_SCREENING_THRESHOLDS


@dataclass
class Input:
    slab_filename: str
    molecule_filename: str
    jobscript_path: str
    submit_command: str

    E_slab_mol : Optional[list[float]]
    jobname_prefix: Optional[str]
    jobscript_ml_path: Optional[str]
    submit_command_ml: Optional[str]



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
    include_surrounding_sites: bool = True
    surrounding_sites_deltaz: float = 1.5
    surrounding_exclude_main: bool = True
    cn_plain_fixed_radius: float = 1.5

@dataclass
class Adsorption_sites:
    mode : str    
    selected_sites: Optional[list[int]]
    high_symmetry_params: Optional[High_symmetry_params]
    coord_number_params: Optional[Coord_number_params]
    surface_height: float = 0.9

@dataclass
class Molecule:
    molecule_axis: dict
    selected_atom_index: int
    x_rot_angles: list[float]
    y_rot_angles: list[float]
    z_rot_angles: list[float]

    individual_rotations: Optional[list[list[float]]]

    vertical_angles: str | bool | list[float] = 'x'
    adsorption_distance_mode: str = 'value'
    target_distance: float = 2.0
    min_distance: float = 1.5
    
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


@dataclass
class Calculation_parameters:
    dft_program: str

    dft_parameters: dict

    def __post_init__(self):
        if self.dft_program.lower() not in SUPPORTED_PROGRAMS:
            raise ValueError('dft_program must be either VASP or QE.')
    

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
        
        self.calculation_parameters = from_dict(data_class=Calculation_parameters, 
                                                data=self.settings_dict["Calculation_parameters"], 
                                                config=Config(type_hooks={str: str.lower}, strict=True))
        
        #read energies of slab and mol if requested
        if READ_ENERGIES and self.input.E_slab_mol is None:
            self.read_E_slab_mol(VERBOSE)

                    



    def read_E_slab_mol(self, VERBOSE : bool = True):
        try:
            from ase.io import read
            slab_en = read(self.input.slab_filename).get_potential_energy()
            mol_en = read(self.input.molecule_filename).get_potential_energy()
            self.input.E_slab_mol = [slab_en, mol_en]
        except:
            self.input.E_slab_mol = [0.0, 0.0]
            if VERBOSE:
                print('It was not possible to obtain slab and molecule energy in any way.' \
                        'Total energies will be shown instead of adsorption energies.')



settings = Settings("../settings.toml")

print(settings.structure.molecule.x_rot_angles)

