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
class CalculationFilesInfo(JsonableBase):
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
        return self.__dict__

    @classmethod
    def fromdict(cls, dct: dict) -> 'CalculationFilesInfo':
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
    bonds : list[BondInfo]

    trajectory : list[AtomsCustom]
    adsorption_energy_evol: list[float]
    final_dz: float

    __xsorb_objtype__ = 'CalculationResults'

    def todict(self) -> dict:
        return self.__dict__

    @classmethod
    def fromdict(cls, dct: dict) -> 'CalculationResults':
        """
        Create an instance of the class from a dictionary.
        Used by xsorb to reconstruct objects after reading from JSON or database.
        """
        # Convert nested objects
        dct['atoms'] = AtomsCustom.fromdict(dct['atoms'])
        dct['trajectory'] = [AtomsCustom.fromdict(atoms) for atoms in dct['trajectory']]
        dct['bonds'] = [BondInfo.fromdict(bond) for bond in dct['bonds']]
        return cls(**dct)


@dataclass
class AdsorptionCalculation(JsonableBase):
    """
    Class that contains all the information about an adsorption calculation.
    """

    adsorption_structure: AdsorptionStructure
    filesinfo: Optional[CalculationFilesInfo]
    results: Optional[CalculationResults]

    __xsorb_objtype__ = 'AdsorptionCalculation'

    def todict(self) -> dict:
        return self.__dict__

    @classmethod
    def fromdict(cls, dct: dict) -> 'AdsorptionCalculation':
        """
        Create an instance of the class from a dictionary.
        Used by xsorb to reconstruct objects after reading from JSON or database.
        """
        # Convert nested objects
        dct['adsorption_structure'] = AdsorptionStructure.fromdict(dct['adsorption_structure'])
        if dct['filesinfo'] is not None:
            dct['filesinfo'] = CalculationFilesInfo.fromdict(dct['filesinfo'])
        if dct['results'] is not None:
            dct['results'] = CalculationResults.fromdict(dct['results'])
        return cls(**dct)
