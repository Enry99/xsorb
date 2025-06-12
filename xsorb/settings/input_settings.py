'''
Module containing all the dataclasses included in the input Settings class,
except for DFTParams, which is defined in a separate module (to keep evrything
related to DFT codes in a single module).
'''

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional

from ase.data import chemical_symbols

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

    @dataclass
    class _RangeSelectionParams:
        '''
        Small helper class to store the settings
        that defines the range selection for the coordination number.
        '''
        mode: str
        value: int | float

        def __post_init__(self):
            if self.mode not in ['max', 'offset']:
                raise ValueError('range_selection mode must be either max or offset.')

    cn_method: str
    atomic_species: Optional[list[str]]
    range_selection: _RangeSelectionParams = field(
        default_factory=lambda: CoordNumberParams._RangeSelectionParams('offset', 2)
    )
    include_surrounding_sites: bool = False
    surrounding_sites_deltaz: float = 1.5
    cn_plain_fixed_radius: float = 1.5

    def __post_init__(self):
        if self.cn_method not in ['plain', 'minimumdistancenn', 'crystalnn']:
            raise ValueError('cn_method must be either plain, MinimumDistanceNN or CystalNN.')
        if self.atomic_species:
            for i, species in enumerate(self.atomic_species):
                self.atomic_species[i] = species.capitalize()
            for species in self.atomic_species:
                if species not in chemical_symbols:
                    raise ValueError(f'atomic_species {species} is invalid.')

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

    @dataclass
    class _MoleculeAxis:
        '''
        Small helper class to store the settings
        that defines the axis of the molecule.
        '''
        mode: str
        values: list[int | float]

        def __post_init__(self):
            if self.mode not in ['atom_indices', 'vector']:
                raise ValueError('molecule_axis mode must be either atom_indices or vector.')
            if self.mode == 'atom_indices' and len(self.values) != 2:
                raise ValueError('molecule_axis values must be a list of two atom indices \
                                    when mode is atom_indices.')
            if self.mode == 'vector' and len(self.values) != 3:
                raise ValueError('molecule_axis values must be a list of three floats \
                                    when mode is vector.')


    molecule_axis: _MoleculeAxis
    selected_atom_indexes: list[int]
    x_rot_angles: list[float]
    y_rot_angles: list[float]
    z_rot_angles: list[float]

    individual_rotations: Optional[list[list[float]]]

    vertical_angles: str | list[float] | None = 'x'
    adsorption_distance_mode: str = 'value'
    target_distance: float = 2.0
    min_distance: float = 1.5
    radius_scale_factor: float = 1.1

    def __post_init__(self):
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
