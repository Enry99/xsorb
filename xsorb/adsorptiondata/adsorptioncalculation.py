'''
Module containing the AdsorptionCalculation class.
'''

from __future__ import annotations
from dataclasses import dataclass, asdict
from typing import Optional, Any

from dacite import from_dict, Config

from xsorb.adsorptiondata.base import JsonableBase, dict_without_none
from xsorb.adsorptiondata.adsorptionstructure import AdsorptionStructure
from xsorb.ase_custom import AtomsCustom

ALLOWED_STATUSES = ('completed', 'incomplete', 'scf_nonconverged')

@dataclass
class CalculationInfo:
    '''
    Small dataclass to store info about the written systems
    '''

    calc_id: str #index(number) or 'slab'/'mol'
    in_file_path: str
    out_file_path: str
    log_file_path: str
    _status: str = 'incomplete' # 'completed', 'incomplete', 'scf_nonconverged'


    __xsorb_objtype__ : str = 'CalculationInfo'

    def __post_init__(self) -> None:
        """
        Post-initialization to ensure that the status is set correctly.
        """

        if self.status not in ALLOWED_STATUSES:
            raise ValueError(f"Status must be one of {ALLOWED_STATUSES}.")

    # make status a property to ensure it is always set correctly
    @property
    def status(self) -> str:
        return self._status

    @status.setter
    def status(self, value: str) -> None:
        if value not in ALLOWED_STATUSES:
            raise ValueError(f"Status must be one of {ALLOWED_STATUSES}.")
        self._status = value


    def db_keys(self) -> dict:
        """
        Returns a dictionary with the keys to be explicitly stored in the database.
        """
        return {
            'calc_id': int(self.calc_id) if self.calc_id.isdigit() else self.calc_id,
            'in_file_path': self.in_file_path,
            'out_file_path': self.out_file_path,
            'log_file_path': self.log_file_path,
            'status': self.status,
        }


@dataclass
class BondInfo:
    """
    Dataclass to store information about a bond between a molecule and a slab.
    """

    mol_atom_id: int
    slab_atom_id: int
    mol_atom_species: str
    slab_atom_species: str
    length: float

    __xsorb_objtype__ : str = 'BondInfo'

    def __str__(self) -> str:
        return f"{self.mol_atom_species}{self.mol_atom_id}-"\
            f"{self.slab_atom_species}{self.slab_atom_id}({self.length:.2f})"


@dataclass
class CalculationResults:
    '''
    Dataclass to store the results of a calculation
    '''

    atoms: AtomsCustom
    adsorption_energy: float
    adsorption_energy_evol: list[float]
    final_dz: float

    bonds : list[BondInfo] | None = None
    trajectory : list[AtomsCustom] | None = None


    __xsorb_objtype__ : str = 'CalculationResults'


    def db_keys(self) -> dict:
        """
        Returns a dictionary with the keys to be explicitly stored in the database.
        """

        dct : dict[str, Any] = {
            'adsorption_energy': self.adsorption_energy,
            'final_dz': self.final_dz,
        }
        if self.bonds is not None:
            dct['bonds'] = ','.join(str(bond) for bond in self.bonds)

        return dct


@dataclass
class AdsorptionCalculation(JsonableBase):
    """
    Class that contains all the information about an adsorption calculation.
    Used to pack all the three main components of the calculation:
    - AdsorptionStructure: the structure of the slab and the molecule
    - CalculationInfo: information about the files used in the calculation
    - CalculationResults: the results of the calculation

    """

    adsorption_structure: AdsorptionStructure
    calc_info: Optional[CalculationInfo] = None
    calc_results: Optional[CalculationResults] = None

    __xsorb_objtype__ : str = 'AdsorptionCalculation'


    def db_keys(self) -> dict:
        """
        Returns a dictionary with the keys to be explicitly stored in the database.
        """

        # merge db_keys from nested objects
        dct = self.adsorption_structure.db_keys()
        if self.calc_info is not None:
            dct.update(self.calc_info.db_keys())
        if self.calc_results is not None:
            dct.update(self.calc_results.db_keys())

        return dct


    def todict(self) -> dict:
        return asdict(self, dict_factory=dict_without_none)

    @classmethod
    def fromdict(cls, dct: dict) -> 'AdsorptionCalculation':
        """
        Create an instance of the class from a dictionary.
        Used by xsorb to reconstruct objects after reading from JSON or database.
        """

        def convert_generators_to_lists(obj):
            if isinstance(obj, dict): #if dict, call this function recursively for each element
                return {k: convert_generators_to_lists(v) for k, v in obj.items()}
            elif isinstance(obj, (list, tuple)): #if list or tuple, call on each element
                return [convert_generators_to_lists(item) for item in obj]
            elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes)):
                # generator, excluding strings and bytes
                try:
                    return list(obj)
                except: #pylint: disable=bare-except
                    return obj
            else:
                return obj

        return from_dict(
            data_class=cls,
            data=convert_generators_to_lists(dct),
            config=Config(type_hooks={AtomsCustom: lambda data: AtomsCustom(data)})
        )
