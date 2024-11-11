#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Enrico Pedretti

'''
Module that contains the class AdsorptionStructuresGenerator,
used to generate the adsorption structures
'''

from __future__ import annotations
from dataclasses import asdict

import numpy as np
from ase import Atoms
from ase.geometry import get_distances
from ase.data import covalent_radii
from ase.data.vdw_alvarez import vdw_radii

from xsorb.ase_custom import AtomsCustom
from xsorb.io.settings import Settings
from xsorb.structures.molecule import Molecule
from xsorb.structures.slab import Slab
from xsorb.structures.properties import (AdsorptionSite, AdsorptionSiteAmorphous,
    MoleculeRotation, AdsorptionStructure, SurroundingSite)

TOL = 0.01 #tolerance for checking if an atom is outside the cell (scaled positions)


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

    def __init__(self, slab: AtomsCustom, mol: AtomsCustom,
                 settings: Settings, verbose : bool = False) -> None:


        self.settings = settings
        self.molecule_rotations = None

        #Slab import from file
        if verbose:
            print('Loading slab...')
        self.slab = Slab(slab=slab,
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
        self.mol = Molecule(mol=mol,
                    atom_indexes=settings.structure.molecule.selected_atom_indexes,
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
                                     which_index : int = 0,
                                     adsite : AdsorptionSiteAmorphous | None = None,
                                     save_image : bool = False,
                                     verbose : bool = True):
        '''
        Generate the molecule rotations, depending on the mode.

        Args:
        - rot_mode: 'standard' for the standard mode, 'surrounding' for the surrounding mode
        - which_index: index of the atom to be used as reference for the rotations in the list of indexes
        - adsite: for the surrounding mode, AdsorptionSiteAmorphous object that contains info to
            rotate the molecule towards the surrounding sites
        - save_image: save an image of the molecule rotations
        - verbose: print additional information during the generation
        '''

        structure_settings = self.settings.structure

        if rot_mode == 'standard':
            if self.molecule_rotations is not None:
                return self.molecule_rotations
            z_rot_angles : list[float | SurroundingSite] = structure_settings.molecule.z_rot_angles
            surrounding_exclude_main = False
        else:
            if adsite is None:
                raise ValueError("adsite must be provided for surrounding mode")
            z_rot_angles = adsite.surrounding_sites
            surrounding_exclude_main = \
                structure_settings.adsorption_sites.coord_number_params.surrounding_exclude_main
            verbose = False

        molecule_rotations = self.mol.generate_molecule_rotations(
            which_index=which_index,
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
        Returns the vertical translation required so that the molecule is not compenetrated
        below the slab, and so that the each mol-slab pair is at least min_distance apart,
        where min_distance can be a value, the sum of the covalent radii, or the sum of
        the van der Waals radii. It works by iteratively translating the molecule upwards
        until the condition is satisfied.

        Args:
        - mol: Atoms object for the molecule
        '''

        mode = self.settings.structure.molecule.adsorption_distance_mode
        slab = self.slab.slab_ase
        min_distance = self.settings.structure.molecule.min_distance
        mult = self.settings.structure.molecule.radius_scale_factor

        if mode == 'value':
            radii = np.array([min_distance]*len(covalent_radii))/2 #dist. is diameter
        elif mode == 'covalent_radius':
            radii = covalent_radii*mult
        elif mode == 'vdw_radius':
            radii = vdw_radii*mult
        else:
            raise ValueError("mode must be one of 'value', 'covalent_radius', 'vdw_radius'")

        dz_tot = 0
        molcopy = mol.copy()

        niter = 0
        while True:

            # First, find the closest slab-molecule pair, not in absolute sense but
            # according to their covalent or van der Waals radii.
            d_matrix = get_distances(slab.positions,molcopy.positions,slab.cell,pbc=True)[1]

            diff_matrix = d_matrix.copy()
            for id1, at1 in enumerate(slab):
                for id2, at2 in enumerate(molcopy):
                    diff_matrix[id1,id2] -= mult*(radii[at1.number] + radii[at2.number])

            i_slab, j_mol, *_ = np.unravel_index(np.argmin(diff_matrix), diff_matrix.shape)

            if diff_matrix[i_slab,j_mol] > -1e-6: #avoid infinite loop due to numerical errors
                break #already above target distance

            #Then, translate the mol upwards so that the z coordinate of the mol atom
            # is equal to the closest slab atom. In this way we can exploit the right triangle
            h = molcopy[j_mol].z - slab[i_slab].z
            if h < 0:
                molcopy.translate([0,0,h])
                dz_tot += h
                h = 0

            #Finally, translate so that the distance between the two atoms is equal
            # to the target distance
            target_dist = mult * (radii[slab[i_slab].number] + radii[molcopy[j_mol].number])
            distance = get_distances(slab[i_slab].position,
                                     molcopy[j_mol].position,
                                     slab.cell,pbc=True)[1][0,0]
            b = np.sqrt(distance**2 - h**2)
            dz = np.sqrt(target_dist**2 - b**2) - h

            molcopy.translate([0,0,dz])
            dz_tot += dz

            niter += 1

            if niter > 100:
                raise RuntimeError("Too many iterations in _mindistance_deltaz")

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
        distance = 0.

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
            atoms : AtomsCustom = mol + self.slab.slab_ase
            mol_indices = list(range(len(mol)))
        else:
            atoms = self.slab.slab_ase + mol
            mol_indices = list(range(len(self.slab.slab_ase), len(self.slab.slab_ase)+len(mol)))
        atoms.cell = self.slab.slab_ase.cell
        atoms.pbc = self.slab.slab_ase.pbc

        return AdsorptionStructure(atoms=atoms,
                                   adsite=adsite,
                                   mol_rot=mol_rot,
                                   distance=distance,
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


    def generate_adsorption_structures(self, write_sites : bool = False,
                                       save_image : bool = False,
                                       verbose : bool = True):
        '''
        Generate all adsorption structures considering the combinations of
        molecular rotations and adsorption sites.

        Args:
        - write_sites: write the adsorption sites to file. This allows to mantain
            the correct site labeling even if ordering is changed after modifying the
            parameters in a following run performed to add sites.
            Use only when launching the calculations, not when just plotting image.
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

        #write sites, so that they are stored for following calculations
        if write_sites:
            Slab.write_sites(adsites)

        #Adsorption of molecule on all adsorption sites for all molecule orientations
        if verbose:
            print('Generating adsorption structures...')

        if sites_settings.mode == 'coord_number'\
            and sites_settings.coord_number_params.include_surrounding_sites:
            rot_mode = 'surrounding'
        else:
            rot_mode = 'standard'

        adsorption_structures : list[AdsorptionStructure] = []
        for mol_ref_index in self.mol.reference_atom_indices:
            for adsite in adsites:
                molecule_rotations = self._generate_molecule_rotations(rot_mode=rot_mode,
                                                                    which_index=mol_ref_index,
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


        #Exclude configs outside cell if requested
        if self.settings.structure.misc.inside_only:
            filtered_adsorption_structures = []
            for ads in adsorption_structures:
                flatposlist = ads.atoms[ads.mol_indices].get_scaled_positions(wrap=False).flatten()
                if np.all(0-TOL < flatposlist) and np.all(flatposlist < 1+TOL):
                    filtered_adsorption_structures.append(ads)
            adsorption_structures = filtered_adsorption_structures

        if verbose:
            print('Adsorption structures generated.')

        return adsorption_structures
