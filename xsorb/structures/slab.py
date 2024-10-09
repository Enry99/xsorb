#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Created on Tue 28 Feb 2023
#@author: Enrico Pedretti

"""
Module that contains the Slab class, with the methods to find adsorption sites

"""

import warnings
from pathlib import Path

import numpy as np
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.coord import find_in_coord_list_pbc
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.local_env import MinimumDistanceNN, CrystalNN
from ase.build.tools import sort
from ase.constraints import FixCartesian
from ase.neighborlist import NeighborList, natural_cutoffs
from ase.geometry.geometry import get_layers

from xsorb.io.utils import ase_custom_read as read
from xsorb.visualize.geometry import save_adsites_image
from xsorb.structures.properties import AdsorptionSite, AdsorptionSiteCrystal, \
    AdsorptionSiteAmorphous, SurroundingSite


class Slab:
    '''
        Class to read slab from file (e.g. Quantum ESPRESSO pwi/pwo or VASP POSCAR),
        find adsorption sites on the surface, and generate the adsorption structures with a molecule
        by placing the molecule on all the different sites

        Initialization parameters:
        - slab_filename: file containing the structure of the slab
        - layers_threshold: deltaz for the atoms to be considered as part of the same layer.
            Used to identify layers when fixing atoms by layer.
        - surface_thickness: thickness (delta_z) of the surface layer (in A)
            to be considered as surface atoms
        - fixed_layers_slab: list of layers of the slab to be fixed
            (counting starts from the bottom, beginning with 0)
        - fixed_indices_slab: list of specific atoms to be fixed (-1: fix all)
            (indices start from 0, with the ordering of the atoms in the input file)
        - fix_slab_xyz: which coordinates to fix for the fixed atoms,
            e.g. [True, True, False] = fix motion in x,y, free to move along z.
        - sort_atoms_by_z: sort the atoms of the slab according to their z coordinate,
            in reverse order (from higher to lower)
        - translate_slab_from_below_cell_bottom: translate slab to make sure that it is
            at least 1 Angstrom from the bottom of the cell
    '''

    def __init__(self, slab_filename : str,
                 surface_thickness : float = 0.9,
                 layers_threshold : float = 0.5,
                 fixed_layers_slab : list | None = None,
                 fixed_indices_slab : list | None = None,
                 fix_slab_xyz : list | None = None,
                 sort_atoms_by_z : bool = False,
                 translate_slab_from_below_cell_bottom : bool = True):

        self.slab_ase = read(filename=slab_filename)
        self.natoms   = len(self.slab_ase)
        self.surface_thickness = surface_thickness

        # First, check if the cell is defined according to the right-handed rule, if not,
        # swap the first two vectors so that the surface normal is found correctly by
        # AdsorbateSiteFinder. Otherwise, the sites will be found on the bottom surface.
        # This also ensures that when using VASP it does not complain.
        cell = self.slab_ase.cell
        if np.dot(np.cross(cell[0], cell[1]), cell[2]) < 0:
            self.slab_ase.set_cell(cell[[1, 0, 2], :])


        if translate_slab_from_below_cell_bottom:
            #translate slab so that the bottom layer is at least 1 angstrom from the bottom
            self.slab_ase.translate([0,0,1-min(self.slab_ase.positions[:,2])])


        #SET CONSTRAINTS (before sorting)##############################
        fixed_atoms_indices = []
        if fixed_layers_slab: # identify the atoms belonging to the various layers
            layer_indicization = get_layers(atoms=self.slab_ase, miller=(0,0,1),
                                            tolerance=layers_threshold)[0]
            fixed_atoms_indices = [i for i, layer in enumerate(layer_indicization) \
                                   if layer in fixed_layers_slab]

        elif fixed_indices_slab:
            if -1 in fixed_indices_slab: #fix all atoms
                fixed_atoms_indices = list(range(self.natoms))
            fixed_atoms_indices = fixed_indices_slab

        c = [FixCartesian(atom_index, mask=[not x for x in fix_slab_xyz]) \
             for atom_index in fixed_atoms_indices]  #negate: in qe 0 = fix, here 1(true)=fix
        self.slab_ase.set_constraint(c) #if no user-defined constraints, c is empty.
        #This is necessary to clean possible constraints read from file
        # (in Molecule there is an if else.)
        ###############################################################

        if sort_atoms_by_z: #sort atoms by height (from higher to lower)
            self.slab_ase = sort(self.slab_ase, tags= -self.slab_ase.positions[:, 2])


        #create pymatgen version of the slab, and initialize the AdsorbateSiteFinder
        with warnings.catch_warnings():
            #to suppress the warning about constraints not supported in pymatgen
            warnings.simplefilter("ignore")
            self.slab_pymat = AseAtomsAdaptor.get_structure(self.slab_ase)

        self.asf = AdsorbateSiteFinder(self.slab_pymat, height=surface_thickness)


    def find_adsorption_sites(self, mode: str, **kwargs) -> \
        list[AdsorptionSiteCrystal | AdsorptionSiteAmorphous]:
        '''
        Find the adsorption sites on the surface of the slab, using the specified mode.

        Args:
        - mode: method to find the adsorption sites. It can be 'high_symmetry' or 'coord_number'.
        - kwargs: additional parameters for the method chosen.

        Returns:
        - list of AdsorptionSite objects, containing the information about the adsorption sites.
        '''

        modes = {'high_symmetry': self.find_adsorption_sites_high_symmetry,
                 'coord_number': self.find_adsorption_sites_coord_number,}

        if mode not in modes:
            raise ValueError(f"mode must be one of {modes.keys()}")

        return modes[mode](**kwargs)


    #TODO: write this file when launching the calculations
    @staticmethod
    def existing_sites() -> list[AdsorptionSite]:
        """
        Read already existing sites from previous calculations, stored in a sites.npy file

        Returns:
        - list of AdsorptionSite objects
        """

        if Path('adsites.npy').is_file():
            existing_sites : list = np.load('adsites.npy', allow_pickle=True).tolist()
            return existing_sites
        else:
            return []


    def _get_symmetrically_equivalent_sets(self, coords_set : list, threshold : float =1e-6):
        """Classifies the adsorption sites into sets of symmetrically equivalent sites.

        Args:
            - coords_set: coordinate set in Cartesian coordinates
            - threshold: tolerance for distance equivalence, used
                as input to in_coord_list_pbc for dupl. checking

        Returns:
            - list of lists of equivalent sites in Cartesian coordinates
        """
        surf_sg = SpacegroupAnalyzer(self.asf.slab, 0.1)
        symm_ops = surf_sg.get_symmetry_operations()
        unique_coords = []
        equivalent_sets = []
        # Convert to fractional
        coords_set = [self.asf.slab.lattice.get_fractional_coords(coords) for coords in coords_set]
        for coords in coords_set:
            in_coord = False
            for op in symm_ops:
                idxs = find_in_coord_list_pbc(unique_coords, op.operate(coords), atol=threshold)
                if idxs:
                    if len(idxs) > 1:
                        raise RuntimeError('ERROR: a new site is symmetrically equivalent' \
                              ' to two INEQUIVALENT sites.')
                    in_coord = True
                    equivalent_sets[idxs[0]] += [self.asf.slab.lattice.get_cartesian_coords(coords)]
                    break
            if not in_coord:
                unique_coords += [coords]
                equivalent_sets += [[self.asf.slab.lattice.get_cartesian_coords(coords)]]

        return equivalent_sets


    def _get_classified_high_symmetry_sites(self, site_types : tuple[str, str, str],
                                           adsites_coords : dict,
                                           selected_sites : list | None = None):
        """
        Build the list of AdsorptionSiteCrystal objects,
        classified according to the site type, and include the info.

        Args:
        - site_types: tuple containing the types of adsorption sites to be classified
        - adsites_coords: dictionary containing the coordinates of the ontop,
            hollow and bridge sites
        - selected_sites (optional): indices of the sites for manual selection

        Returns:
        - list of AdsorptionSiteCrystal objects,
        containing the information about the adsorption sites.
        """

        #create structure containing only surface atoms
        surf = self.asf.slab.copy().remove_sites(self.asf.slab.subsurface_sites())
        #flatten the surface (z=0 for all)
        for i in range(len(surf)): surf[i].z = 0 # pylint: disable=consider-using-enumerate,multiple-statements

        nn = MinimumDistanceNN(tol=0.2) #increased tol to identify as 3-fold the sites that are
        #at the center of a non-perfeclty equilater triangle


        #create a list with all the sites, with info for each one
        all_adsites : list[AdsorptionSiteCrystal] = []

        #handle the case of existing sites
        existing_sites = self.existing_sites()
        all_adsites.extend(existing_sites)

        i_site = len(all_adsites)
        for site_type, sites_coords_list in zip(site_types,
            [adsites_coords['ontop'], adsites_coords['hollow'], adsites_coords['bridge']]):

            for site_coords in sites_coords_list:

                unique_id = "{0:.2f},{1:.2f}".format(*site_coords[:2]) #pylint: disable=consider-using-f-string

                if selected_sites and i_site not in selected_sites:
                    continue
                if any(unique_id == x.unique_id for x in existing_sites):
                    continue


                #dummy structure just to place one atom in the site and obtain CN and NN list
                dummy_coords = site_coords.copy()
                dummy_coords[2] = 0.2 #place the dummy atom just above the surface level z=0
                surf.append('O', dummy_coords, coords_are_cartesian=True)
                coord_n = nn.get_cn(surf, len(surf)-1)
                nn_list = nn.get_nn(surf, len(surf)-1)
                surf.remove_sites([len(surf)-1]) #remove dummy atom


                # add further information to the site types
                if site_type == 'ontop':
                    first_nn_species = nn_list[0].site.species_string
                    info = f"{site_type} {first_nn_species}"
                elif site_type == 'hollow':
                    info = f"{site_type} {coord_n}-fold"
                elif site_type == 'bridge':
                    if coord_n>=4: #attemps to fix the problem of fake bridges for 4-fold sites
                        info = f"hollow {coord_n}-fold"
                    else:
                        if len(nn_list) >=2:
                            distance = np.linalg.norm(nn_list[0].coords[:2] - nn_list[1].coords[:2])
                            info = f"{site_type} {distance}"
                        else:
                            info = f"{site_type}"
                else:
                    raise ValueError('Invalid site type.')

                all_adsites.append(AdsorptionSiteCrystal(label=i_site,
                                                         coords=site_coords,
                                                         info=info))
                i_site += 1

        return all_adsites


    def find_adsorption_sites_high_symmetry(self, symm_reduce : float = 0.01,
                              no_obtuse_hollow : bool = True,
                              selected_sites : list | None = None,
                              save_image : bool = False,
                              figname : str = 'adsorption_sites.png',
                              verbose : bool = False):
        '''
        Find adsorption sites using high symmetry mode.
        Optionally save a figure with the sites on the surface.

        Args:
            -symm_reduce: Pymatgen's AdsorbateSiteFinder.find_adsorption_sites parameter.
                It is a theshold for removing symmetrically equivalent sites.
            -no_obtuse_hollow: Pymatgen's AdsorbateSiteFinder.find_adsorption_sites parameter.
                Avoid considering hollow sites inside obtuse triangles of the Delaunay triangulation
                of topmost layer used to find sites.
            -selected_sites: indices of the sites to be returned by this function,
                selected between those found by AdsorbateSiteFinder
            -save_image: decide wether to save a png image of the sites
            -figname: filename of the image.
        '''

        if verbose:
            print('Finding adsorption sites...')

        site_types = ('ontop', 'hollow', 'bridge')

        #adsites is a dict which contains the ontop, hollow and bridge sites
        adsites_coords = self.asf.find_adsorption_sites(distance=0,
                                                 symm_reduce=0,
                                                 no_obtuse_hollow=no_obtuse_hollow).pop('all')


        #find the equivalent sets of sites according to symm_reduce.
        #in each set, selec the site closest to te cell center as the representative site
        for site_type in site_types:
            #sort the sites by distance from the center of the cell
            sorted_sites = sorted(adsites_coords[site_type],
                key=lambda x: np.linalg.norm(
                    x - self.asf.slab.lattice.get_cartesian_coords([0.5, 0.5, 0]))
                    )

            #find the equivalent sets of sites
            equivalent_sets_list = self._get_symmetrically_equivalent_sets(sorted_sites,symm_reduce)

            #select the representative site for each set and assign the sites to the adsites dict
            adsites_coords[site_type]=[equivalent_set[0] for equivalent_set in equivalent_sets_list]

        #classify the adsorption sites, building the list of AdsorptionSiteCrystal objects
        all_adsites = self._get_classified_high_symmetry_sites(site_types,
                                                               adsites_coords,
                                                               selected_sites)


        if verbose:
            print('Adsorption sites found.')

        if save_image: #save png to visualize the identified sites
            save_adsites_image(mode='high_symmetry',
                               adsites=all_adsites,
                               slab_pymat=self.slab_pymat,
                               figname=figname,
                               verbose=verbose)

        return all_adsites


    def _coord_number_surface_analysis(self,
                                      cn_method: str,
                                      cn_plain_fixed_radius: float | None = 1.5,
                                      verbose: bool=False):
        """
        Analyzes the coordination number of surface sites in a slab structure.

        Args:

        - cn_method : str
            The method to use for coordination number calculation. Options are:
            'plain' - Uses a fixed radius or natural cutoffs.
            'minimumdistancenn' - Uses pymatgen's MinimumDistanceNN.
            'crystalnn' - Uses pymatgen's CrystalNN with weights.
        - cn_plain_fixed_radius : float or None, optional
            Fixed radius for the 'plain' method. If set to None, natural cutoffs are used.
            Default is 1.5.
        - verbose : bool, optional
            If True, prints additional information during the calculation. Default is False.

        Returns:
        - surf_coords : list
            List of coordinates of the surface sites.
        - cn_list : list
            List of coordination numbers for the surface sites.
        - surf_sites_indices : list
            List of indices of the surface sites in the slab structure.

        """

        #Select which sites belong to the surface. Sites whose z coordinate
        # is below self.surface_thickness from the maximum z of the atoms within
        # a circle of radius 2.5 A around the atom are considered surface sites.
        # This allows to have a more accurate selection of the surface sites
        # in case of non-flat surfaces.
        circle_radius = 2.5
        surf_coords = []
        surf_sites_indices = []
        for atom in self.slab_ase:
            #find the max z of atoms within a circle of radius 2.5 A around the atom
            max_z = np.max([atom.position[2] for atom in self.slab_ase \
                    if np.linalg.norm(atom.position[:2] - atom.position[:2])  <circle_radius])
            if atom.position[2] > max_z - self.surface_thickness:
                surf_coords.append(atom.position)
                surf_sites_indices.append(atom.index)


        if cn_method == 'plain':
            #standard coordination number (integer no. of nearest neighbours)
            #with fixed radius. If radius is set to None, the natural_cutoffs are used.
            if verbose:
                print('Using ase.neighborlist {0} to find coordination numbers.'.format( # pylint: disable=consider-using-f-string
                    'with fixed radius for all atoms' if cn_plain_fixed_radius \
                          else 'with ase.neighborlist.natural_cutoffs.'))


            if cn_plain_fixed_radius:
                cutoffs = [cn_plain_fixed_radius]*len(self.slab_ase)
            else:
                cutoffs = natural_cutoffs(self.slab_ase, mult=1.1)
            nl = NeighborList(cutoffs, skin=0, self_interaction=False, bothways=True)
            nl.update(self.slab_ase)
            cm = nl.get_connectivity_matrix()
            cn_list = [cm[idx].count_nonzero() for idx in surf_sites_indices]

        elif cn_method == 'minimumdistancenn':
            if verbose:
                print('Using pymatgen.analysis.local_env.MinimumDistanceNN' \
                        ' to find coordination numbers.')

            nn = MinimumDistanceNN(tol=0.2)
            cn_list = [nn.get_cn(self.asf.slab, idx) for idx in surf_sites_indices]

        elif cn_method == 'crystalnn':
            if verbose:
                print('Using pymatgen.analysis.local_env.CrystalNN with weigths' \
                        ' to find coordination numbers.')


            nn = CrystalNN(weighted_cn=True)
            with warnings.catch_warnings(): #to suppress the warning
                warnings.simplefilter("ignore")
                cn_list = [nn.get_cn(self.asf.slab, idx, use_weights=True) \
                           for idx in surf_sites_indices]

        else:
            raise ValueError('Invalid cn_method.')

        return surf_coords, cn_list, surf_sites_indices


    def _add_surrounding_adsites(self, adsites : list[AdsorptionSiteAmorphous],
                                surf_sites_indices : list[int],
                                surrounding_sites_deltaz: float = 1.5):
        """
        Add (inplace) the surrounding sites to the adsites list, for each adsite.

        Args:
        - adsites: list of AdsorptionSiteAmorphous objects
        - surf_sites_indices: indices of the surface sites in the slab
        - surrounding_sites_deltaz: deltaz from the main site to consider a site as surrounding

        """

        nn = MinimumDistanceNN(tol=0.2)


        all_main_idxs = [adsite.unique_id for adsite in adsites]
        all_connected_idxs = []

        for adsite in adsites:

            adsite.surrounding_sites = []

            nnsites = nn.get_nn_info(self.asf.slab, adsite.unique_id)

            nn_counter_for_this_site = 0
            for nnsite in nnsites:

                nnindex = nnsite['site_index']

                if (nnindex in surf_sites_indices if surrounding_sites_deltaz is None \
                    else abs(adsite.coords[2]-nnsite['site'].coords[2]) < surrounding_sites_deltaz):

                    atom_species = self.asf.slab[nnindex].species_string

                    #calculate vector from main site to surrounding site
                    vector = nnsite['site'].coords - adsite.coords

                    surrounding_site = SurroundingSite(
                        label = f'{adsite.unique_id}.{nn_counter_for_this_site}',
                        coords = nnsite['site'].coords,
                        info = atom_species,
                        atom_index = nnindex,
                        duplicate_surrounding = nnindex in all_connected_idxs,
                        duplicate_main = nnindex in all_main_idxs,
                        vector = vector
                    )
                    adsite.surrounding_sites.append(surrounding_site)
                    all_connected_idxs.append(nnindex)
                    nn_counter_for_this_site += 1


    def _get_classified_coord_number_sites(self,
                                           surf_coords : list,
                                          cn_list : list,
                                          surf_sites_indices : list,
                                          max_cn : float,
                                          selected_sites : list | None = None,
                                          selected_atomic_species : list | None = None,
                                          include_surrounding_sites : bool = False,
                                          surrounding_sites_deltaz: float | None = 1.5):
        """
        Build the list of AdsorptionSiteAmorphous objects, and include the info.
        """

        #create a list with all the sites, with info for each one
        all_adsites : list[AdsorptionSite] = []

        #handle the case of existing sites
        existing_sites = self.existing_sites()
        all_adsites.extend(existing_sites)


        i_site = len(all_adsites)
        for coords, cn, idx in zip(surf_coords, cn_list, surf_sites_indices):

            #here unique id is the atom index in the slab, since we are considering atoms as sites
            unique_id = idx
            atom_species = self.asf.slab[idx].species_string

            if selected_sites and i_site not in selected_sites:
                continue
            if selected_atomic_species and atom_species not in selected_atomic_species:
                continue
            if any(unique_id == x.unique_id for x in existing_sites):
                continue

            if cn <= max_cn:
                info = f"{atom_species}(cn={cn:.1f})"
                all_adsites.append(AdsorptionSiteAmorphous(label=str(i_site),
                                                           coords=coords,
                                                           info=info,
                                                           atom_index=idx,
                                                           coordination_number=cn))
                i_site += 1

        if include_surrounding_sites:
            self._add_surrounding_adsites(
                adsites=all_adsites,
                surf_sites_indices=surf_sites_indices,
                surrounding_sites_deltaz=surrounding_sites_deltaz)

        return all_adsites


    def find_adsorption_sites_coord_number(self,
                                        cn_method: str,
                                        cn_plain_fixed_radius : float | None = 1.5,
                                        max_cn_offset : float | None = 2,
                                        max_cn: float | None = None,
                                        selected_sites: list | None = None,
                                        selected_atomic_species: list | None = None,
                                        include_surrounding_sites : bool = False,
                                        surrounding_sites_deltaz: float | None = 1.5,
                                        save_image: bool = False,
                                        figname: str = 'adsorption_sites.png',
                                        verbose: bool=False):
        '''
        Find adsorption sites using coordination number mode.

        Args:
        - cn_method: method to calculate the coordination number. It can be
            'plain', 'minimumdistancenn' or 'crystalnn'.
        - cn_plain_fixed_radius: fixed radius for the coordination number calculation,
            if cn_method is 'plain'.
        - max_cn_offset: the max. coord. numb. will be set to min(cn) + max_cn_offset.
        - max_cn: maximum coordination number allowed. Priority over max_cn_offset.
        - selected_sites: indices of the sites to be returned by this function,
            selected between those found by AdsorbateSiteFinder.
        - selected_atomic_species: list of atomic species to be considered as adsorption sites.
        - include_surrounding_sites: decide wether to include the surrounding sites in the output.
        - surrounding_sites_deltaz: deltaz from the main site to consider a site as surrounding.
        - save_image: decide wether to save a png image of the sites.
        - figname: filename of the image.
        -verbose: decide wether to print messages during the calculation.

        Returns:
        - list of AdsorptionSite objects, containing the information about the adsorption sites.
        '''

        if verbose:
            print('Finding adsorption sites...')

        surf_coords, cn_list, surf_sites_indices = \
            self._coord_number_surface_analysis(cn_method, cn_plain_fixed_radius, verbose)

        max_cn = max_cn if max_cn is not None else min(cn_list) + max_cn_offset

        #build the list of AdsorptionSiteAmorphous objects, and include the info
        all_adsites = self._get_classified_coord_number_sites(
            surf_coords=surf_coords,
            cn_list=cn_list,
            surf_sites_indices=surf_sites_indices,
            max_cn=max_cn,
            selected_sites=selected_sites,
            selected_atomic_species=selected_atomic_species,
            include_surrounding_sites=include_surrounding_sites,
            surrounding_sites_deltaz=surrounding_sites_deltaz)

        if verbose:
            print('Adsorption sites found.')


        if save_image: #save png to visualize the identified sites
            save_adsites_image(
                mode='coord_number_surrounding' if include_surrounding_sites else 'coord_number',
                adsites=all_adsites,
                slab_pymat=self.slab_pymat,
                figname=figname,
                verbose=verbose)

        return all_adsites
