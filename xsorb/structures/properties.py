#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Enrico Pedretti

'''
Module to store all the dataclasses that describe adsorption structures, sites and rotations

'''

from __future__ import annotations
from dataclasses import dataclass

from xsorb.ase_custom.atoms import AtomsCustom


@dataclass
class MoleculeRotation:
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
    - reference_atom: int, index of the reference atom in the moelcule

    Properties:
    - unique_id: string that fully identifies the rotation

    Can be compared with the equality operator, that compares the unique_id.
    '''
    atoms: AtomsCustom
    xrot: str
    yrot: str
    zrot: str
    mol_atom: int

    @classmethod
    def from_dict(cls, data: dict):
        '''
        Creates a MoleculeRotation object from a dictionary
        '''
        return cls(AtomsCustom(data['atoms']),
                   data['xrot'],
                   data['yrot'],
                   data['zrot'],
                   data['mol_atom'])

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


@dataclass
class AdsorptionSite:
    '''
    Base class of Adsorption Site, to be inherited by the two different modes.

    Contains:
    - label: str, numeric label of the site as it appears in the adsorption sites figure.
    - coords: list[float], x,y,z coordinates of the site
    - info: str, additional information about the site

    Properties:
    - unique_id: string that fully identifies the site

    Can be compared with the equality operator, that compares the unique_id.
    '''

    label: str
    coords: list[float]
    info: str

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

    #define equality as the equality of the unique_id
    def __eq__(self, other) -> bool:
        if not isinstance(other, AdsorptionSite):
            return NotImplemented
        return self.unique_id == other.unique_id


@dataclass
class AdsorptionSiteCrystal(AdsorptionSite):
    '''
    Class to store the information of an adsorption site in the high-symmetry mode.

    Contains:
    - label: str, numeric label of the site as it appears in the adsorption sites figure.
    - coords: list[float], x,y,z coordinates of the site
    - info: str, additional information about the site, e.g. "ontop Cu" or "hollow 3-fold"

    Properties:
    - unique_id: string that fully identifies the site

    Can be compared with the equality operator, that compares the unique_id.
    '''

    type: str #type of the site, e.g. ontop, bridge, hollow, etc.

    #for now leave it empty, we might need to add more attributes in the future


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


    def __str__(self) -> str:
        return self.label


@dataclass
class AdsorptionStructure:
    '''
    Class to store the information of an adsorption structure

    Contains:
    - atoms: Atoms object of the adsorption structure
    - adsite: AdsorptionSite object of the adsorption site
    - mol_rot: MoleculeRotation object of the rotated molecule
    - distance: float, distance between the reference atom of the molecule
    - mol_indices: list[int], indices of the atoms of the molecule
    and the adsorption site

    Methods:
    - to_info_dict: returns a dictionary with the information of the AdsorptionStructure object
    '''
    atoms: AtomsCustom
    adsite: AdsorptionSite
    mol_rot: MoleculeRotation
    distance : float
    mol_indices: list[int]


    @classmethod
    def from_dict(cls, data: dict):
        '''
        Creates an AdsorptionStructure object from a dictionary
        '''
        return cls(AtomsCustom(data['atoms']),
                   AdsorptionSite.from_dict(data['adsite']),
                   MoleculeRotation.from_dict(data['mol_rot']),
                   data['distance'],
                   data['mol_indices'])


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
        return ("site", "site_info", "mol_atom", "initial_dz", "xrot", "yrot", "zrot")

    def to_info_dict(self):
        """
        Returns a dictionary with the information of the AdsorptionStructure object
        """

        try: #to make the ase db check happy
            site = int(self.adsite.label)
        except ValueError:
            site = self.adsite.label
        try: #to make the ase db check happy
            zrot = float(self.mol_rot.zrot)
        except ValueError:
            zrot = self.mol_rot.zrot
        infodict = {"site": site,
                    "site_info": self.adsite.info,
                    "mol_atom": self.mol_rot.mol_atom,
                    "initial_dz": self.distance,
                    "xrot": float(self.mol_rot.xrot),
                    "yrot": float(self.mol_rot.yrot),
                    "zrot": zrot}

        #check that the keys are the same as the column names.
        #(to avoid introducing bugs when adding new columns in the code)
        assert self.dataframe_column_names() == tuple(infodict.keys())

        return infodict
