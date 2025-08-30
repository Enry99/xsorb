#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Enrico Pedretti

'''
Module to store all the dataclasses that describe adsorption structures, sites and rotations.
Classes shoul all inherit from JsonableBase to be able to be serialized to JSON.

'''

from __future__ import annotations
from dataclasses import dataclass

import numpy as np

from xsorb.ase_custom.atoms import AtomsCustom
from xsorb.adsorptiondata.base import JsonableBase


@dataclass
class MoleculeRotation(JsonableBase):
    '''
    Class to store a rotated molecule and the rotation angles.
    The rotation angles are stored as strings, to be able to store also the
    SurroundingSite object when using the coordination number method.
    Angles are not intended to reproduce the rotation of the molecule given
    the initial structure, since the rotated molecule is already stored in
    the Atoms object.

    Contains:
    - atoms: Atoms object of the rotated molecule
    - xrot: string with the x rotation angle
    - yrot: string with the y rotation angle
    - zrot: string with the z rotation angle
    - mol_atom: int, index of the reference atom in the moelcule

    Properties:
    - unique_id: string that fully identifies the rotation

    Can be compared with the equality operator, that compares the unique_id.
    '''

    atoms: AtomsCustom
    xrot: str
    yrot: str
    zrot: str
    mol_atom: int

    __xsorb_objtype__ : str = "MoleculeRotation"


    def db_keys(self) -> dict:
        """
        Returns a dictionary with the keys to be explicitly stored in the database.
        """

        def is_numeric(s):
            try:
                float(s)
                return True
            except (ValueError, TypeError):
                return False

        xrot = self.xrot if not is_numeric(self.xrot) else float(self.xrot)
        yrot = self.yrot if not is_numeric(self.xrot) else float(self.yrot)
        zrot = self.zrot if not is_numeric(self.zrot) else float(self.zrot)
        return {
            'xrot': xrot,
            'yrot': yrot,
            'zrot': zrot,
            'mol_atom': self.mol_atom
        }

    @property
    def unique_id(self):
        '''
        String that fully identifies the rotation
        '''
        return f"{self.xrot},{self.yrot},{self.zrot}"

    #define equality as the equality of the unique_id
    def __eq__(self, other) -> bool:
        if not isinstance(other, MoleculeRotation):
            return NotImplemented
        return self.unique_id == other.unique_id

    def todict(self):
        """Convert class instance to a dictionary."""
        return self.__dict__

    @classmethod
    def fromdict(cls, dct: dict):
        '''
        Creates a MoleculeRotation object from a dictionary
        '''
        dct['atoms'] = AtomsCustom(dct['atoms'])
        return cls(**dct)


@dataclass
class AdsorptionSite(JsonableBase):
    '''
    Base class of Adsorption Site, to be inherited by the two different modes.

    Contains:
    - label: str, label of the site as it appears in the adsorption sites figure, e.g. "1" or "2.1"
    - coords: list[float], x,y,z coordinates of the site
    - info: str, additional information about the site

    Properties:
    - unique_id: string that fully identifies the site

    Can be compared with the equality operator, that compares the unique_id.
    '''

    label: str
    coords: list[float]
    info: str


    def db_keys(self) -> dict:
        """
        Returns a dictionary with the keys to be explicitly stored in the database.
        """
        return {
            'site': self.label if not self.label.isdigit() else int(self.label),
            #'coords': self.unique_id,
            'site_info': self.info
        }

    @classmethod
    def from_dict(cls, data: dict):
        '''
        Creates an AdsorptionSite object from a dictionary
        '''
        return cls(data['label'], list(data['coords']), data['info'])

    @property
    def unique_id(self):
        '''
        String that fully identifies the site
        '''
        return "{0:.2f},{1:.2f},{2:.2f}".format(*self.coords) #pylint: disable=consider-using-f-string

    def __eq__(self, other) -> bool:
        if not isinstance(other, AdsorptionSite):
            return NotImplemented
        return np.allclose(self.coords, other.coords)

    def todict(self):
        """Convert class instance to a dictionary."""
        return self.__dict__

    @classmethod
    def fromdict(cls, dct: dict):
        '''
        Creates an AdsorptionSite object from a dictionary
        '''
        return cls(**dct)


@dataclass
class AdsorptionSiteCrystal(AdsorptionSite):
    '''
    Class to store the information of an adsorption site in the high-symmetry mode.

    Contains:
    - label: str, numeric label of the site as it appears in the adsorption sites figure.
    - coords: list[float], x,y,z coordinates of the site
    - info: str, additional information about the site, e.g. "ontop Cu" or "hollow 3-fold"
    - type: str, type of the site, e.g. "ontop", "bridge", "hollow", etc.

    Properties:
    - unique_id: string that fully identifies the site

    Can be compared with the equality operator, that compares the unique_id.
    '''

    type: str #type of the site, e.g. ontop, bridge, hollow, etc.

    __xsorb_objtype__ : str = "AdsorptionSiteCrystal"

    def todict(self):
        return self.__dict__

    @classmethod
    def fromdict(cls, dct: dict):
        '''
        Creates an AdsorptionSiteCrystal object from a dictionary
        '''
        return cls(**dct)

    def __eq__(self, other) -> bool:
        if not isinstance(other, AdsorptionSiteCrystal):
            return NotImplemented
        return super().__eq__(other)


@dataclass
class AdsorptionSiteAmorphous(AdsorptionSite):
    '''
    Class to store the information of an adsorption site in the amorphous mode.

    Contains:
    - label: str, numeric label of the site as it appears in the adsorption sites figure.
    - coords: list[float], x,y,z coordinates of the site
    - info: str, additional information about the site, e.g. Cu(cn=3)
    - atom_index: int, index of the atom in the Atoms object
    - coordination_number: float, coordination number of the site
    - surrounding_sites: list[SurroundingSite], list of the surrounding sites

    Properties:
    - unique_id: string that fully identifies the site

    '''

    atom_index: int
    coordination_number: float | None = None
    surrounding_sites: list['SurroundingSite'] | None = None

    __xsorb_objtype__ : str = "AdsorptionSiteAmorphous"

    def todict(self):
        #remove None values from the dictionary to write into ase db
        dct = self.__dict__.copy()
        dct = {k: v for k, v in dct.items() if v is not None}
        return dct

    @classmethod
    def fromdict(cls, dct: dict):
        '''
        Creates an AdsorptionSiteAmorphous object from a dictionary
        '''
        if 'surrounding_sites' in dct:
            dct['surrounding_sites'] = [SurroundingSite.fromdict(ss) \
                    for ss in dct['surrounding_sites']]
        return cls(**dct)

    def __eq__(self, other) -> bool:
        if not isinstance(other, AdsorptionSiteAmorphous):
            return NotImplemented
        return super().__eq__(other)


@dataclass
class SurroundingSite(AdsorptionSite):
    '''
    Class to store the information of a surrounding site in the amorphous mode.

    Contains:
    - label: str, numeric label of the site as it appears in the adsorption sites figure.
    - coords: list[float], x,y,z coordinates of the site
    - info: str, additional information about the site, e.g. Cu
    - atom_index: int, index of the atom in the Atoms object
    - duplicate_surrounding: bool, if the site is a duplicate of the surrounding sites
    - duplicate_main: bool, if the site is a duplicate of the main sites
    - vector: list[float], vector from the main site to the surrounding site

    Properties:
    - unique_id: string that fully identifies the site

    '''

    atom_index: int
    duplicate_surrounding: bool
    duplicate_main: bool
    vector: list[float]   #vector from the main site to the surrounding site

    __xsorb_objtype__ : str = "SurroundingSite"

    def __str__(self) -> str:
        # to be printed in csvfile as z_rot
        return f"to_{self.label}"

    def todict(self):
        return self.__dict__

    @classmethod
    def fromdict(cls, dct: dict):
        '''
        Creates a SurroundingSite object from a dictionary
        '''
        return cls(**dct)


@dataclass
class AdsorptionStructure(JsonableBase):
    '''
    Class to store the information of an adsorption structure

    Contains:
    - atoms: Atoms object of the adsorption structure
    - adsite: AdsorptionSite object of the adsorption site
    - mol_rot: MoleculeRotation object of the rotated molecule
    - distance: float, distance between the reference atom of the molecule and the adsorption site
    - mol_indices: list[int], indices of the atoms of the molecule

    Methods:
    - to_info_dict: returns a dictionary with the information of the AdsorptionStructure object
    '''

    atoms: AtomsCustom
    adsite: AdsorptionSite
    mol_rot: MoleculeRotation
    distance : float
    mol_indices: list[int]

    __xsorb_objtype__ : str = "AdsorptionStructure"


    def db_keys(self) -> dict:
        """
        Returns a dictionary with the keys to be explicitly stored in the database.
        """
        dct = self.adsite.db_keys()
        dct.update(self.mol_rot.db_keys())
        dct.update({'initial_dz': self.distance})

        #check that the keys are the same as the column names.
        #(to avoid introducing bugs when adding new columns in the code)
        assert set(self.dataframe_column_names()) == set(dct.keys())

        return dct

    @property
    def slab_indices(self):
        '''
        Returns the indices of the atoms of the slab
        '''
        return [i for i in range(len(self.atoms)) if i not in self.mol_indices]


    @staticmethod
    def dataframe_column_names():
        '''
        Returns the names of the columns of the AdsorptionStructure object
        '''
        return ("site", "site_info", "xrot", "yrot", "zrot", "mol_atom", "initial_dz")


    def todict(self):
        return self.__dict__


    @classmethod
    def fromdict(cls, dct: dict):
        '''
        Creates an AdsorptionStructure object from a dictionary
        '''

        dct['atoms'] = AtomsCustom(dct['atoms'])
        dct['mol_rot'] = MoleculeRotation.fromdict(dct['mol_rot'])

        #handle different types of AdsorptionSite
        if dct['adsite']['__xsorb_objtype__'] == "AdsorptionSiteCrystal":
            dct['adsite'] = AdsorptionSiteCrystal.fromdict(dct['adsite'])
        elif dct['adsite']['__xsorb_objtype__'] == "AdsorptionSiteAmorphous":
            dct['adsite'] = AdsorptionSiteAmorphous.fromdict(dct['adsite'])
        elif dct['adsite']['__xsorb_objtype__'] == "SurroundingSite":
            dct['adsite'] = SurroundingSite.fromdict(dct['adsite'])
        else:
            raise ValueError(f"Unknown AdsorptionSite type: {dct['adsite']['__xsorb_objtype__']}")

        return cls(**dct)
