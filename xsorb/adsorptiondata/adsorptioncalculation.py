'''
Module containing the AdsorptionCalculation class.
'''

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional

from xsorb.adsorptiondata.base import JsonableBase
from xsorb.adsorptiondata.adsorptionstructure import AdsorptionStructure
from xsorb.ase_custom import AtomsCustom


@dataclass
class CalculationInfo(JsonableBase):
    '''
    Small dataclass to store info about the written systems
    '''

    calc_id: int | str #index or 'slab'/'mol'
    in_file_path: str
    out_file_path: str
    log_file_path: str
    job_id: int | None = None
    job_status : str | None = None # 'running', 'completed', 'failed', 'cancelled', None

    __xsorb_objtype__ = 'CalculationFilesInfo'

    def todict(self) -> dict:
        dct = self.__dict__.copy()
        dct = {k: v for k, v in dct.items() if v is not None}
        return dct

    @classmethod
    def fromdict(cls, dct: dict) -> 'CalculationInfo':
        """
        Create an instance of the class from a dictionary.
        Used by xsorb to reconstruct objects after reading from JSON or database.
        """
        return cls(**dct)


@dataclass
class BondInfo(JsonableBase):
    """
    Dataclass to store information about a bond between a molecule and a slab.
    """

    mol_atom_id: int
    slab_atom_id: int
    mol_atom_species: str
    slab_atom_species: str
    length: float

    __xsorb_objtype__ = 'BondInfo'

    def __str__(self) -> str:
        return f"{self.mol_atom_species}{self.mol_atom_id}-"\
            f"{self.slab_atom_species}{self.slab_atom_id}({self.length:.2f})"

    def todict(self) -> dict:
        return self.__dict__

    @classmethod
    def fromdict(cls, dct: dict) -> 'BondInfo':
        """
        Create an instance of the class from a dictionary.
        Used by xsorb to reconstruct objects after reading from JSON or database.
        """
        return cls(**dct)


@dataclass
class CalculationResults(JsonableBase):
    '''
    Dataclass to store the results of a calculation
    '''

    atoms: AtomsCustom
    adsorption_energy: float
    status : str #'completed', 'incomplete'
    scf_nonconverged : bool
    adsorption_energy_evol: list[float]
    final_dz: float

    bonds : list[BondInfo] | None
    trajectory : list[AtomsCustom] | None


    __xsorb_objtype__ = 'CalculationResults'

    def todict(self) -> dict:
        dct = self.__dict__.copy()
        dct = {k: v for k, v in dct.items() if v is not None}
        return dct

    @classmethod
    def fromdict(cls, dct: dict) -> 'CalculationResults':
        """
        Create an instance of the class from a dictionary.
        Used by xsorb to reconstruct objects after reading from JSON or database.
        """
        # Convert nested objects
        dct['atoms'] = AtomsCustom.fromdict(dct['atoms'])
        if 'trajectory' in dct:
            dct['trajectory'] = [AtomsCustom.fromdict(atoms) for atoms in dct['trajectory']]
        if 'bonds' in dct:
            dct['bonds'] = [BondInfo.fromdict(bond) for bond in dct['bonds']]
        return cls(**dct)


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
    calc_info: Optional[CalculationInfo]
    calc_results: Optional[CalculationResults]

    __xsorb_objtype__ = 'AdsorptionCalculation'

    def todict(self) -> dict:
        dct = self.__dict__.copy()
        dct = {k: v for k, v in dct.items() if v is not None}
        return dct

    @classmethod
    def fromdict(cls, dct: dict) -> 'AdsorptionCalculation':
        """
        Create an instance of the class from a dictionary.
        Used by xsorb to reconstruct objects after reading from JSON or database.
        """
        # Convert nested objects
        dct['adsorption_structure'] = AdsorptionStructure.fromdict(dct['adsorption_structure'])
        if dct['calc_info'] is not None:
            dct['calc_info'] = CalculationInfo.fromdict(dct['calc_info'])
        if dct['calc_results'] is not None:
            dct['calc_results'] = CalculationResults.fromdict(dct['calc_results'])
        return cls(**dct)
