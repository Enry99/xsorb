'''
Module that contains the class AdsorptionStructuresGenerator,
used to generate the adsorption structures
'''

from dataclasses import asdict

import numpy as np
from ase import Atoms
from ase.geometry import get_distances
from ase.data import covalent_radii
from ase.data.vdw_alvarez import vdw_radii

from xsorb.settings import Settings
from xsorb.structures.molecule import Molecule
from xsorb.structures.slab import Slab
from xsorb.structures.properties import AdsorptionSite, AdsorptionSiteAmorphous, \
    MoleculeRotation, AdsorptionStructure



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
                    surface_thickness=settings.structure.adsorption_sites.surface_thickness,
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


    def _mindistance_deltaz(self, mol: Atoms):
        '''
        Returns the vertical translation required so that the closest mol-slab atom pair is at least
        at least min_distance away from each other along the z axis, where min_distance can be
        a value, the sum of the covalent radii, or the sum of the van der Waals radii.
        It works by iteratively translating the molecule upwards until the condition is satisfied.

        Args:
        - mol: Atoms object for the molecule
        '''

        mode = self.settings.structure.molecule.adsorption_distance_mode
        slab = self.slab.slab_ase
        min_distance = self.settings.structure.molecule.min_distance
        mult = self.settings.structure.molecule.radius_scale_factor

        dz_tot = 0
        molcopy = mol.copy()

        while True:
            #Iteratively translate the molecule upwards until the condition is satisfied
            #The iterative procedure is necessary if part of the molecule is below the slab

            vector_matrix, d_matrix = get_distances(slab.positions,mol.positions,slab.cell,pbc=True)
            candidate_translations = []
            for id1, at1 in enumerate(slab):
                for id2, at2 in enumerate(mol):

                    distance = d_matrix[id1, id2]
                    h = np.dot(vector_matrix[id1, id2], [0,0,1])
                    b = np.sqrt(distance**2 - h**2)

                    if mode == 'covalent_radius':
                        D = mult*(covalent_radii[at1.number] + covalent_radii[at2.number]) #pylint: disable=invalid-name
                    elif mode == 'vdw_radius':
                        D = mult*(vdw_radii[at1.number] + vdw_radii[at2.number]) #pylint: disable=invalid-name
                    else:
                        D = min_distance #pylint: disable=invalid-name

                    if distance < D:
                        dz = np.sqrt(D**2 - b**2) - h
                        candidate_translations.append(dz)

            if candidate_translations:
                dz = np.max(candidate_translations)
                dz_tot += dz
                molcopy.translate([0,0,dz])
            else:
                break

        return dz_tot


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
        distance = 0

        #place the molecule in the adsorption site, then translate it upwards by the target height
        mol.translate(adsite.coords)

        if self.settings.structure.molecule.adsorption_distance_mode == 'value':
            #Use the target distance from the settings. Otherwise just let the iterative
            #min_z_distance_from_surf function handle the distance
            distance += self.settings.structure.molecule.target_distance
            mol.translate([0,0,distance])

        #Check for min_z_distance_from_surf and translate accordingly
        further_transl = self._mindistance_deltaz(mol)
        if further_transl:
            mol.translate( [0, 0, further_transl] )
            distance += further_transl


        if self.settings.structure.misc.mol_before_slab:
            atoms : Atoms = mol + self.slab.slab_ase
            slab_indices = list(range(len(mol), len(mol)+len(self.slab.slab_ase)))
            mol_indices = list(range(len(mol)))
        else:
            atoms : Atoms = self.slab.slab_ase + mol
            slab_indices = list(range(len(self.slab.slab_ase)))
            mol_indices = list(range(len(self.slab.slab_ase), len(self.slab.slab_ase)+len(mol)))
        atoms.cell = self.slab.slab_ase.cell

        return AdsorptionStructure(atoms=atoms,
                                   adsite=adsite,
                                   mol_rot=mol_rot,
                                   distance=distance,
                                   slab_indices=slab_indices,
                                   mol_indices=mol_indices)


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
