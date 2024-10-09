#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue 28 Feb 2023
@author: Enrico Pedretti

Small helper class to handle the slab

"""

import warnings
from pathlib import Path

import numpy as np
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.coord import find_in_coord_list_pbc
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.local_env import MinimumDistanceNN, CrystalNN
from ase.io import read
from ase.build.tools import sort
from ase.constraints import FixCartesian
from ase.neighborlist import NeighborList, natural_cutoffs
from ase.geometry.geometry import get_layers

from xsorb.visualize.geometry import save_adsites_image
from xsorb.structures.properties import AdsorptionSite, SurroundingSite
from xsorb import ase_custom








class Slab:
    '''
        Class to read slab from file (e.g. Quantum ESPRESSO pwi/pwo or VASP POSCAR),
        find adsorption sites on the surface, and generate the adsorption structures with a molecule
        by placing the molecule on all the different sites

        Initialization parameters:
        - slab_filename: file containing the structure of the slab
        - layers_threshold: deltaz for the atoms to be considered as part of the same layer. 
        Used to identify layers when fixing atoms by layer.
        - surface_sites_height: deltaz from the topmost atom to identify which atoms are surface atoms,
        used by Pymatgen to find adsorption sites (it is the 'height' parameter of Pymatgen's AdsorbateSiteFinder)
        - fixed_layers_slab: list of layers of the slab to be fixed 
        (counting starts from the bottom, beginning with 0)
        - fixed_indices_slab: list of specific atoms to be fixed (-1: fix all)
        (indices start from 0, with the ordering of the atoms in the input file)
        - fix_slab_xyz: which coordinates to fix for the fixed atoms, e.g. [True, True, False] = fix motion in x,y, free to move along z.
        - sort_atoms_by_z: sort the atoms of the slab according to their z coordinate, in reverse order (from higher to lower)
        - translate_slab_from_below_cell_bottom: translate slab to make sure that it is at least 1 Angstrom from the bottom of the cell
    '''

    def __init__(self, slab_filename : str, 
                 surface_sites_height : float = 0.9, 
                 layers_threshold : float = 0.5, 
                 fixed_layers_slab : list | None = None, 
                 fixed_indices_slab : list | None = None, 
                 fix_slab_xyz : list | None = None, 
                 sort_atoms_by_z : bool = False,
                 translate_slab_from_below_cell_bottom : bool = True):
        
        self.slab_ase = ase_custom.Atoms_custom(read(filename=slab_filename))
        self.natoms   = len(self.slab_ase)

        # First, check if the cell is defined according to the right-handed rule, if not,
        # swap the first two vectors so that the surface normal is found correctly by AdsorbateSiteFinder. 
        # Otherwise, the sites will be found on the bottom surface.
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
            layer_indicization = get_layers(atoms=self.slab_ase, miller=(0,0,1), tol=layers_threshold)[0]
            fixed_atoms_indices = [i for i, layer in enumerate(layer_indicization) if layer in fixed_layers_slab]

        elif fixed_indices_slab:
            if -1 in fixed_indices_slab: #fix all atoms
                fixed_atoms_indices = list(range(self.natoms))
            fixed_atoms_indices = fixed_indices_slab

        c = [FixCartesian(atom_index, mask=[not x for x in fix_slab_xyz]) \
             for atom_index in fixed_atoms_indices]  #we need to negate: in qe 0 = fix, here 1(true)=fix
        self.slab_ase.set_constraint(c) #if no user-defined constraints, c is empty. 
        #This is necessary to clean possible constraints read from file (in Molecule there is an if else.)
        ###############################################################

        if sort_atoms_by_z: #sort atoms by height (from higher to lower)
            self.slab_ase = sort(self.slab_ase, tags= -self.slab_ase.positions[:, 2])  

        
        #create pymatgen version of the slab, and initialize the AdsorbateSiteFinder
        with warnings.catch_warnings(): #to suppress the warning about constraints not supported in pymatgen
            warnings.simplefilter("ignore")
            self.slab_pymat = AseAtomsAdaptor.get_structure(self.slab_ase)

        self.asf = AdsorbateSiteFinder(self.slab_pymat, height=surface_sites_height) 


    def find_adsorption_sites(self, mode: str, **kwargs) -> list[AdsorptionSite]:   

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
        """

        if Path('adsites.npy').is_file():
            existing_sites : list = np.load('adsites.npy', allow_pickle=True).tolist()
            return existing_sites
        else:
            return []


    def get_symmetrically_equivalent_sets(self, coords_set : list, threshold : float =1e-6):
        """Classifies the adsorption sites into sets of symmetrically equivalent sites.
    
        Args:
            coords_set: coordinate set in Cartesian coordinates
            threshold: tolerance for distance equivalence, used
                as input to in_coord_list_pbc for dupl. checking
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
                        print('ERROR: a new site is symmetrically equivalent to two INEQUIVALENT sites.')
                    in_coord = True
                    equivalent_sets[idxs[0]] += [self.asf.slab.lattice.get_cartesian_coords(coords)]
                    break
            if not in_coord:
                unique_coords += [coords]
                equivalent_sets += [[self.asf.slab.lattice.get_cartesian_coords(coords)]]

        return equivalent_sets
    

    def find_adsorption_sites_high_symmetry(self, symm_reduce : float = 0.01,  
                              no_obtuse_hollow : bool = True, 
                              selected_sites : list | None = None,   
                              SAVE_IMAGE : bool = False,
                              figname : str = 'adsorption_sites.png',
                              VERBOSE : bool = False):
        '''
        Returns a list of cartesian coordinates of the adsites, and a list of labels ('ontop', x, y).
        Optionally it saves a figure with the sites on the surface.

        Args:
            -symm_reduce: Pymatgen's AdsorbateSiteFinder.find_adsorption_sites parameter. It is a theshold for removing symmetrically equivalent sites.
            -no_obtuse_hollow: Pymatgen's AdsorbateSiteFinder.find_adsorption_sites parameter. Avoid considering hollow sites inside obtuse triangles of the Delaunay triangulation of topmost layer used to find sites.
            -selected_sites: indices of the sites to be returned by this function, selected between those found by AdsorbateSiteFinder
            -save_image: decide wether to save a png image of the sites
            -figname: filename of the image.
        '''

        if VERBOSE: print('Finding adsorption sites...')

        site_types = ('ontop', 'hollow', 'bridge')

        #adsites is a dict which contains the ontop, hollow and bridge sites
        adsites = self.asf.find_adsorption_sites(distance=0, 
                                                 symm_reduce=0, 
                                                 no_obtuse_hollow=no_obtuse_hollow).pop('all')
        
        
        #find the equivalent sets of sites according to symm_reduce.
        #in each set, selec the site closest to te cell center as the representative site
        for site_type in site_types:
            #sort the sites by distance from the center of the cell
            sorted_sites = sorted(adsites[site_type], 
                key=lambda x: np.linalg.norm(x - self.asf.slab.lattice.get_cartesian_coords([0.5, 0.5, 0])))
            
            #find the equivalent sets of sites
            equivalent_sets_list = self.get_symmetrically_equivalent_sets(sorted_sites, symm_reduce)

            #select the representative site for each set and assign the selected sites to the adsites dict
            adsites[site_type] = [equivalent_set[0] for equivalent_set in equivalent_sets_list]
            
        
        
        #create structure containing only surface atoms
        surf = self.asf.slab.copy().remove_sites(self.asf.slab.subsurface_sites()) #remove subsurface atoms
        for i in range(len(surf)): surf[i].z = 0 #flatten the surface (z=0 for all)

        nn = MinimumDistanceNN(tol=0.2) #increased tol to identify as 3-fold the sites that are at the 
        #center of a non-perfeclty equilater triangle

        
        #create a list with all the sites, with info for each one
        all_adsites : list[AdsorptionSite] = []
        
        #handle the case of existing sites
        existing_sites = self.existing_sites()
        all_adsites.extend(existing_sites)

        i_site = len(all_adsites)
        for site_type, sites_list in zip(site_types, \
                                         [adsites['ontop'], adsites['hollow'], adsites['bridge']]):
            for site in sites_list:

                unique_id = "{0:.2f},{1:.2f}".format(*site[:2])

                if selected_sites and i_site not in selected_sites: continue
                if any(unique_id == x.unique_id for x in existing_sites): continue

                
                #dummy structure just to place one atom in the site and obtain CN and NN list
                coords = site.tolist()
                coords[2] = 0.2 #place the dummy atom just above the surface level z=0
                surf.append('O', coords, coords_are_cartesian=True)
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
                    if(coord_n>=4): #attemps to fix the problem of fake bridges for 4-fold sites
                        info = f"hollow {coord_n}-fold"
                    else:
                        if len(nn_list) >=2:
                            distance = np.linalg.norm(nn_list[0].coords[:2] - nn_list[1].coords[:2])
                            info = f"{site_type} {distance}"
                        else: 
                            info = f"{site_type}"
                else:
                    raise ValueError('Invalid site type.')

                all_adsites.append(AdsorptionSite(label=i_site, 
                                                  coords=site, 
                                                  info=info, 
                                                  unique_id=unique_id))
                i_site += 1
                

        if VERBOSE: print('Adsorption sites found.')

        if(SAVE_IMAGE): #save png to visualize the identified sites
            save_adsites_image(mode='high_symmetry',
                               adsites=all_adsites,
                               slab_pymat=self.slab_pymat, 
                               figname=figname, 
                               VERBOSE=VERBOSE)

        return all_adsites


    def coord_number_surface_analysis(self,
                                   cn_method: str, 
                                   cn_plain_fixed_radius : float | None = 1.5,
                                   VERBOSE: bool=False):
                                   
        surf_coords = [s.coords for s in self.asf.surface_sites]
        surf_sites_indices = [i for i in range(len(self.asf.slab.sites)) \
                              if np.any(np.all(self.asf.slab.cart_coords[i] == surf_coords, axis=1))]

        if cn_method == 'plain':
            #standard coordination number (integer no. of nearest neighbours)
            #with fixed radius. If radius is set to None, the natural_cutoffs are used.
            if VERBOSE: print('Using ase.neighborlist {0} to find coordination numbers.'\
                  .format('with fixed radius for all atoms' if cn_plain_fixed_radius \
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
            if VERBOSE: print('Using pymatgen.analysis.local_env.MinimumDistanceNN \
                              to find coordination numbers.')
            
            nn = MinimumDistanceNN(tol=0.2)
            cn_list = [nn.get_cn(self.asf.slab, idx) for idx in surf_sites_indices]
        
        elif cn_method == 'crystalnn':
            if VERBOSE: print('Using pymatgen.analysis.local_env.CrystalNN with weigths \
                              to find coordination numbers.')
            

            nn = CrystalNN(weighted_cn=True)
            with warnings.catch_warnings(): #to suppress the warning
                warnings.simplefilter("ignore")
                cn_list = [nn.get_cn(self.asf.slab, idx, use_weights=True) for idx in surf_sites_indices]

        else:
            raise ValueError('Invalid cn_method.')
        
        return surf_coords, cn_list, surf_sites_indices


    def add_surrounding_adsites(self, adsites : list[AdsorptionSite], 
                                surf_sites_indices : list[int],
                                surrounding_sites_deltaz: float = 1.5):
        """
        Add (inplace) the surrounding sites to the adsites list, for each adsite.
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
                    else abs(adsite.coords[2] - nnsite['site'].coords[2]) < surrounding_sites_deltaz):

                    atom_species = self.asf.slab[nnindex].species_string

                    #calculate vector from main site to surrounding site
                    vector = nnsite['site'].coords - adsite.coords

                    surrounding_site = SurroundingSite(
                        label = f'{adsite.unique_id}.{nn_counter_for_this_site} {atom_species}',
                        coords = nnsite['site'].coords,
                        duplicate_surrounding = nnindex in all_connected_idxs,
                        duplicate_main = nnindex in all_main_idxs,
                        vector = vector,
                        unique_id = nnindex    
                    )
                    adsite.surrounding_sites.append(surrounding_site)
                    all_connected_idxs.append(nnindex)
                    nn_counter_for_this_site += 1


    def find_adsorption_sites_coord_number(self,               
                                        cn_method: str, 
                                        cn_plain_fixed_radius : float | None = 1.5, 
                                        max_cn_offset : float | None = 2, 
                                        max_cn: float | None = None, 
                                        selected_sites: list | None = None, 
                                        selected_atomic_species: list | None = None,
                                        include_surrounding_sites : bool = False, 
                                        surrounding_sites_deltaz: float | None = 1.5, 
                                        SAVE_IMAGE: bool = False,
                                        figname: str = 'adsorption_sites.png',
                                        VERBOSE: bool=False):

        if VERBOSE: print('Finding adsorption sites...')

        surf_coords, cn_list, surf_sites_indices = \
            self.coord_number_surface_analysis(cn_method, cn_plain_fixed_radius, VERBOSE)
        
        max_cn = max_cn if max_cn is not None else min(cn_list) + max_cn_offset

        
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

            if selected_sites and i_site not in selected_sites: continue
            if selected_atomic_species and atom_species not in selected_atomic_species: continue
            if any(unique_id == x.unique_id for x in existing_sites): continue

            if cn <= max_cn:
                info = f"{atom_species}(cn={cn:.1f})"
                all_adsites.append(AdsorptionSite(label=i_site, 
                                                  coords=coords, 
                                                  type='ontop', 
                                                  info=info, 
                                                  unique_id=unique_id,
                                                  coordination_number=cn))
                i_site += 1

        if include_surrounding_sites:
            self.add_surrounding_adsites(
                adsites=all_adsites, 
                surf_sites_indices=surf_sites_indices,
                surrounding_sites_deltaz=surrounding_sites_deltaz)
            

        if VERBOSE: print('Adsorption sites found.')


        if SAVE_IMAGE: #save png to visualize the identified sites
            save_adsites_image(
                mode='coord_number_surrounding' if include_surrounding_sites else 'coord_number',
                adsites=all_adsites, 
                slab_pymat=self.slab_pymat,
                figname=figname, 
                VERBOSE=VERBOSE)

        return all_adsites

