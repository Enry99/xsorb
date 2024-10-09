'''
Module that contains the class AdsorptionStructuresGenerator, 
used to generate the adsorption structures
'''

from dataclasses import asdict

from ase import Atoms
from ase.data import atomic_numbers, covalent_radii

from xsorb.structures.utils import closest_pair
from xsorb.settings import Settings
from xsorb.structures.molecule import Molecule
from xsorb.structures.slab import Slab
from xsorb.structures.properties import AdsorptionSite, AdsorptionSiteAmorphous, \
    MoleculeRotation, AdsorptionStructure


def mindistance_deltaz(slab : Atoms, mol: Atoms, min_z_distance_from_surf : float):
    '''
    Returns the vertical translation required so that the closest mol-slab atom pair is at least
    min_z_distance_from_surf apart (or half the sum of the covalent radii) along z

    Args:
    - slab: Atoms object for the slab
    - mol: Atoms object for the molecule
    - min_distance: minimum required distance along z
    '''

    dz_tot = 0
    molcopy = mol.copy()
    while True:
        #First, find the closest slab-mol atoms pair
        i_mol, j_slab, _ = closest_pair(slab, molcopy)

        #Find the z coordinates of the closest atoms pair, and half their covalent distance
        half_covalent_distance = 0.5 * (covalent_radii[atomic_numbers[slab[j_slab].symbol]] \
                                    + covalent_radii[atomic_numbers[molcopy[i_mol].symbol]])
        zmol = molcopy[i_mol].position[2]
        zslab = slab[j_slab].position[2]

        #Calculate the distance required to enforce the minimum distance
        necessary_min_z_dist = max(min_z_distance_from_surf, half_covalent_distance)
        if zmol < zslab + necessary_min_z_dist:
            dz = zslab + necessary_min_z_dist - zmol
            dz_tot += dz
            molcopy.translate([0,0,dz])
        else:
            break

    return dz_tot



class AdsorptionStructuresGenerator:
    '''
    Class to generate adsorption structures for a given slab and molecule, 
    obtained as the combinations of the molecule on different adsorption sites
    and different molecular rotations.
    When instantiated, it reads slab and molecule from file, and does the 
    preprocessing of the structural parameters.
    Then the generate_adsorption_structures method can be used to get the 
    adsorption structures.

    Initialization parameters:
    - settings: Settings object containing the input parameters
    - verbose: print additional information during the initialization

    Attributes:
    - settings: Settings object containing the input parameters
    - slab: Slab object containing the slab structure
    - mol: Molecule object containing the molecule structure
    - molecule_rotations: list of MoleculeRotation objects containing the rotated molecules.
        It is initialized as None, and is stored for subsequent calls to avoid recomputing it.

    Methods:
    - generate_adsorption_structures: generate all adsorption structures considering 
        the combinations of molecular rotations and adsorption sites
    '''

    def __init__(self, settings: Settings, verbose : bool = False) -> None:


        self.settings = settings
        self.molecule_rotations = None

        #Slab import from file
        if verbose:
            print('Loading slab...')
        self.slab = Slab(slab_filename=settings.input.slab_filename,
                    surface_sites_height=settings.structure.adsorption_sites.surface_height, 
                    layers_threshold=settings.structure.constraints.layers_height, 
                    fixed_layers_slab=settings.structure.constraints.fixed_layers_slab, 
                    fixed_indices_slab=settings.structure.constraints.fixed_indices_slab, 
                    fix_slab_xyz=settings.structure.constraints.fix_slab_xyz,
                    sort_atoms_by_z=settings.structure.misc.sort_atoms_by_z,
                    translate_slab_from_below_cell_bottom=settings.structure.misc.translate_slab)
        if verbose:
            print('Slab loaded.')


        #Molecule import from file
        if verbose:
            print('Loading molecule...')
        self.mol = Molecule(molecule_filename=settings.input.molecule_filename,
                    atom_index=settings.structure.molecule.selected_atom_index,
                    molecule_axis_atoms=settings.structure.molecule.molecule_axis["values"] \
                        if settings.structure.molecule.molecule_axis["mode"] == "atom_indices" \
                            else None,
                    axis_vector=settings.structure.molecule.molecule_axis["values"] \
                        if settings.structure.molecule.molecule_axis["mode"] == "vector" \
                            else None,
                    fixed_indices_mol=settings.structure.constraints.fixed_indices_mol,
                    fix_mol_xyz=settings.structure.constraints.fix_mol_xyz)
        if verbose:
            print('Molecule loaded.')



    def _generate_molecule_rotations(self,
                                     rot_mode : str = 'standard',
                                     adsite : AdsorptionSiteAmorphous | None = None,
                                     save_image : bool = False,
                                     verbose : bool = True):
        '''
        Generate the molecule rotations, depending on the mode.

        Args:
        - rot_mode: 'standard' for the standard mode, 'surrounding' for the surrounding mode
        - adsite: for the surrounding mode, AdsorptionSiteAmorphous object that contains info to 
            rotate the molecule towards the surrounding sites
        - save_image: save an image of the molecule rotations
        - verbose: print additional information during the generation
        '''

        structure_settings = self.settings.structure

        if rot_mode == 'standard':
            if self.molecule_rotations is not None:
                return self.molecule_rotations
            z_rot_angles = structure_settings.molecule.z_rot_angles
            surrounding_exclude_main = False
        else:
            if adsite is None:
                raise ValueError("adsite must be provided for surrounding mode")
            z_rot_angles = adsite.surrounding_sites
            surrounding_exclude_main = \
                structure_settings.adsorption_sites.coord_number_params.surrounding_exclude_main
            verbose = False

        molecule_rotations = self.mol.generate_molecule_rotations(
            x_rot_angles=structure_settings.molecule.x_rot_angles,
            y_rot_angles=structure_settings.molecule.y_rot_angles,
            z_rot_angles=z_rot_angles,
            vert_angles_list=structure_settings.molecule.vertical_angles,
            individual_rotations=structure_settings.molecule.individual_rotations,
            surrounding_exclude_main=surrounding_exclude_main,
            save_image=save_image,
            verbose=verbose)

        if rot_mode == 'standard':
            #in case of standard, store the rotations to avoid calling the function multiple times
            self.molecule_rotations = molecule_rotations

        return molecule_rotations


    def _put_together_slab_and_mol(self, adsite : AdsorptionSite,
                                  mol_rot : MoleculeRotation):
        '''
        Place the molecule on the adsorption site and translate it upwards by the target height.

        Args:
        - adsite: AdsorptionSite object containing the adsorption site
        - mol_rot: MoleculeRotation object containing the rotated molecule

        Returns:
        - AdsorptionStructure object containing the Atoms object and info 
            on the adsorption site and the molecule rotation
        '''


        mol = mol_rot.atoms.copy()

        #place the molecule in the adsorption site, then translate it upwards by the target height
        mol.translate(adsite.coords)

        distance = self.settings.structure.molecule.target_distance
        mol.translate([0,0,distance])

        #Check for min_z_distance_from_surf and translate accordingly
        further_transl = mindistance_deltaz(self.slab.slab_ase,
                                            mol,
                                            self.settings.structure.molecule.min_distance)
        if further_transl:
            mol.translate( [0, 0, further_transl] )
            distance += further_transl

        atoms : Atoms = mol + self.slab.slab_ase if self.settings.structure.misc.mol_before_slab \
            else self.slab.slab_ase + mol
        atoms.cell = self.slab.slab_ase.cell

        return AdsorptionStructure(atoms=atoms,
                                   adsite=adsite,
                                   mol_rot=mol_rot,
                                   distance=distance)


    def _get_structures_for_vertical_surrounding_sites(self, adsite : AdsorptionSiteAmorphous):
        '''
        Obtain adsorption structures on surrounding sites of an adsorption site in the amorphous
        mode, when the molecule is placed vertically.

        Args:
        - adsite: AdsorptionSiteAmorphous object containing the surrounding sites

        Returns:
        - adsorption_structures: list of AdsorptionStructure objects
        '''

        adsorption_structures : list[AdsorptionStructure] = []

        for surr_site in adsite.surrounding_sites:

            if surr_site.duplicate_main or surr_site.duplicate_surrounding:
                continue

            molecule_rotations = self._generate_molecule_rotations(rot_mode='standard')
            for mol_rot in molecule_rotations:
                #if y_angle != +/- 90, continue
                if int(float(mol_rot.yrot)) not in [90, -90]:
                    continue
                structure = self._put_together_slab_and_mol(adsite=surr_site,
                                                           mol_rot=mol_rot)
                adsorption_structures.append(structure)


        return adsorption_structures


    def generate_adsorption_structures(self, save_image : bool = False, verbose : bool = True):
        '''
        Generate all adsorption structures considering the combinations of 
        molecular rotations and adsorption sites.

        Args:
        - save_image: save an image of the adsorption sites and of the molecular rotations.
        - verbose: print additional information during the generation

        Returns:
        - adsorption_structures: list of AdsorptionStructure objects, each containing 
            the Atoms object and info on the adsorption site and the molecule rotation.
        '''

        sites_settings = self.settings.structure.adsorption_sites


        #Find adsorption sites and labels (site type and x,y coords.)
        mode_params = {
            'high_symmetry': sites_settings.high_symmetry_params,
            'coord_number': sites_settings.coord_number_params,}
        mode = sites_settings.mode
        if mode not in mode_params:
            raise ValueError(f"mode must be one of {mode_params.keys()}")


        adsites = self.slab.find_adsorption_sites(
            mode=mode,
            **asdict(mode_params[mode]),
            selected_sites=sites_settings.selected_sites,
            save_image=save_image,
            verbose=verbose)


        #Adsorption of molecule on all adsorption sites for all molecule orientations
        if verbose:
            print('Generating adsorption structures...')

        if sites_settings.mode == 'coord_number'\
            and sites_settings.coord_number_params.include_surrounding_sites:
            rot_mode = 'surrounding'
        else:
            rot_mode = 'standard'

        adsorption_structures : list[AdsorptionStructure] = []
        for adsite in adsites:
            molecule_rotations = self._generate_molecule_rotations(rot_mode=rot_mode,
                                                                  adsite=adsite,
                                                                  save_image=save_image,
                                                                  verbose=verbose)
            for mol_rot in molecule_rotations:
                structure = self._put_together_slab_and_mol(adsite=adsite,mol_rot=mol_rot)
                adsorption_structures.append(structure)


            #handle the case of vertical molecule on surrounding sites
            if rot_mode == 'surrounding':
                adsorption_structures.extend(
                    self._get_structures_for_vertical_surrounding_sites(adsite))


        if verbose:
            print('Adsorption structures generated.')

        return adsorption_structures
