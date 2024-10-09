'''
Module to store all the dataclasses that describe adsorption structures, sites and rotations

'''

from dataclasses import dataclass

from ase import Atoms


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

    Properties:
    - unique_id: string that fully identifies the rotation

    Can be compared with the equality operator, that compares the unique_id.
    '''
    atoms: Atoms
    xrot: str
    yrot: str
    zrot: str

    @property
    def unique_id(self):
        '''
        String that fully identifies the rotation
        '''
        return f"{self.xrot},{self.yrot},{self.zrot}"

    #define equality as the equality of the unique_id
    def __eq__(self, other : 'MoleculeRotation') -> bool:
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

    @property
    def unique_id(self):
        '''
        String that fully identifies the site
        '''
        return "{0:.2f},{1:.2f},{2:.2f}".format(*self.coords)

    #define equality as the equality of the unique_id
    def __eq__(self, other : 'AdsorptionSite') -> bool:
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
    and the adsorption site

    Properties:
    - dataframe_column_names: tuple with the names of the columns of the AdsorptionStructure object

    Methods:
    - to_dataframe_row: returns a list with the data of the AdsorptionStructure object
    to be used as a row in a pandas DataFrame
    '''
    atoms: Atoms
    adsite: AdsorptionSite
    mol_rot: MoleculeRotation
    distance : float

    @property
    def dataframe_column_names(self):
        """
        Returns a tuple with the names of the columns of the AdsorptionStructure object,
        to be used as column names in a pandas DataFrame
        """
        return ("site",
                "site_info",
                "x",
                "y",
                "z",
                "distance",
                "xrot",
                "yrot",
                "zrot")

    def to_dataframe_row(self):
        """
        Returns a list with the data of the AdsorptionStructure object, 
        to be used as a row in a pandas DataFrame
        """
        return [self.adsite.label,
                self.adsite.info,
                *self.adsite.coords,
                self.distance,
                self.mol_rot.xrot,
                self.mol_rot.yrot,
                self.mol_rot.zrot]

        #df = pd.DataFrame(columns=['a', 'b'])
        #df_line = [1,2]
        #df.loc[0] = df_line
