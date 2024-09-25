#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue 28 Feb 2023
@author: Enrico Pedretti

Small helper class to handle the slab

"""

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patheffects as PathEffects
from pymatgen.analysis.adsorption import AdsorbateSiteFinder, plot_slab, get_rot
from pymatgen.analysis.local_env import MinimumDistanceNN
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core import Structure
from ase import Atoms
from ase.io import read
from ase.build.tools import sort
from ase.data import atomic_numbers, covalent_radii
from ase.neighborlist import NeighborList, natural_cutoffs
from ase.constraints import FixCartesian
from ase.geometry import get_distances
from ase.geometry.geometry import get_layers
from xsorb.visualize.geometry import save_adsites_image
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
                 layers_threshold : float = 0.5, 
                 surface_sites_height : float = 0.9, 
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
            layer_indicization = get_layers(atoms=self.slab_ase, miller=(0,0,1)[0], tol=layers_threshold)
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
        import warnings
        with warnings.catch_warnings(): #to suppress the warning about constraints not supported in pymatgen
            warnings.simplefilter("ignore")
            self.slab_pymat = AseAtomsAdaptor.get_structure(self.slab_ase)

        self.asf = AdsorbateSiteFinder(self.slab_pymat, height=surface_sites_height) 


    def find_adsorption_sites(self, crystal : bool = True, **kwargs):   
        if crystal:
            return self.find_adsorption_sites_crystal(**kwargs)
        else:
            return self.find_adsorption_sites_amorphous(**kwargs)


    def find_adsorption_sites_crystal(self, symm_reduce : float = 0.01, 
                              near_reduce : float = 0.01, 
                              no_obtuse_hollow : bool = True, 
                              selected_sites : list = None,   
                              SAVE_IMAGE : bool = False,
                              figname : str = 'adsorption_sites.png',
                              VERBOSE : bool = False):
        '''
        Returns a list of cartesian coordinates of the adsites, and a list of labels ('ontop', x, y).
        Optionally it saves a figure with the sites on the surface.

        Args:
            -symm_reduce: Pymatgen's AdsorbateSiteFinder.find_adsorption_sites parameter. It is a theshold for removing symmetrically equivalent sites.
            -near_reduce: Pymatgen's AdsorbateSiteFinder.find_adsorption_sites parameter. Threshold for removing sites duplicates (increase it to reduce duplicates).
            -no_obtuse_hollow: Pymatgen's AdsorbateSiteFinder.find_adsorption_sites parameter. Avoid considering hollow sites inside obtuse triangles of the Delaunay triangulation of topmost layer used to find sites.
            -selected_sites: indices of the sites to be returned by this function, selected between those found by AdsorbateSiteFinder
            -save_image: decide wether to save a png image of the sites
            -figname: filename of the image.
        '''

        if VERBOSE: print('Finding adsorption sites...')

        adsites = self.asf.find_adsorption_sites(distance=0, 
                                                 symm_reduce=symm_reduce, 
                                                 near_reduce=near_reduce, 
                                                 no_obtuse_hollow=no_obtuse_hollow)
        if selected_sites:
            sel_adsites = [adsites['all'][i] for i in selected_sites]
        else:
            sel_adsites = adsites['all']

        adsite_labels = []

        #create structure containing only surface atoms
        surf_coords = [s.coords for s in self.asf.surface_sites]
        nonsurf_sites_indices = [i for i in range(len(self.asf.slab.sites)) if not np.any(np.all(self.asf.slab.cart_coords[i] == surf_coords, axis=1))]
        slab = self.asf.slab.copy()
        slab.remove_sites(nonsurf_sites_indices)
        for i in range(len(slab)): slab[i].z = 0 #flatten the surface (z=0 for all)

        nn = MinimumDistanceNN(tol=0.2) #increased tol to identify as 3-fold the sites that are at the center of a non-perfeclty equilater triangle

        #run over all the slab_adsites, classifiying them by checking if the
        #i-th element of 'all' is in one of the three lists 'ontop', 'hollow', 'bridge'        
        for i, site in enumerate(sel_adsites):
            #dummy structure just to place one atom in the site
            coords = site.tolist()
            coords[2] = 0.2 #place the dummy atom just above the surface level z=0
            slab.append('O', coords, coords_are_cartesian=True)
            coord_n = nn.get_cn(slab, len(slab)-1)
            nn_list = nn.get_nn(slab, len(slab)-1) 
            slab.remove_sites([len(slab)-1]) #remove dummy atom     

            # add further information to the site types
            if any((site == x).all() for x in adsites['ontop']):
                first_nn_species = nn_list[0].species_string #atomic species of the atom just below the ontop site
                adsite_labels.append('{0} ontop_{1},{2:.3f},{3:.3f},'.format(i, first_nn_species, *site[:2]))   
            elif any((site == x).all() for x in adsites['hollow']): #coordination number of the hollow site
                adsite_labels.append('{0} hollow_c{1},{2:.3f},{3:.3f},'.format(i, coord_n, *site[:2])) 
            else: 
                if(coord_n>=4): #attemps to fix the problem of fake bridges for 4-fold sites
                    adsite_labels.append('{0} hollow_c{1},{2:.3f},{3:.3f},'.format(i, coord_n, *site[:2]))
                else:
                    if len(nn_list) >=2:
                        distance = np.linalg.norm(nn_list[0].coords[:2] - nn_list[1].coords[:2])
                        adsite_labels.append('{0} bridge_{1:.2f},{2:.3f},{3:.3f},'.format(i, distance, *site[:2]))
                    else: adsite_labels.append('{0} bridge,{1:.3f},{2:.3f},'.format(i, *site[:2]))
        
        
        if VERBOSE: print('Adsorption sites found.')

        if(SAVE_IMAGE): #save png to visualize the identified sites
            save_adsites_image(adsites=sel_adsites, 
                               adsite_labels=adsite_labels, 
                               slab_pymat=self.slab_pymat, 
                               figname=figname, 
                               VERBOSE=VERBOSE)

        return sel_adsites, adsite_labels, None


    def find_adsorption_sites_amorphous(self, **kwargs):

        if kwargs.get('VERBOSE'): print('Finding adsorption sites...')

        surf_coords = [s.coords for s in self.asf.surface_sites]
        surf_sites_indices = [i for i in range(len(self.asf.slab.sites)) if np.any(np.all(self.asf.slab.cart_coords[i] == surf_coords, axis=1))]

        if kwargs.get('cn_method') == 'CrystalNN':
            print('Using CrystalNN to find coordination numbers.')
            from pymatgen.analysis.local_env import CrystalNN
            nn = CrystalNN(weighted_cn=True)
            import warnings
            with warnings.catch_warnings(): #to suppress the warning
                warnings.simplefilter("ignore")
                cn_list = [nn.get_cn(self.asf.slab, idx, use_weights=True) for idx in surf_sites_indices]
        
        elif kwargs.get('cn_method') == 'MinimumDistanceNN':
            print('Using MinimumDistanceNN to find coordination numbers.')
            nn = MinimumDistanceNN(tol=0.2)
            cn_list = [nn.get_cn(self.asf.slab, idx) for idx in surf_sites_indices]
        
        else: #standard coordination number (integer no. of nearest neighbours) with ase.neighborlist
            
            print('Using ase.neighborlist {0}to find coordination numbers.'\
                  .format('with fixed radius for all atoms ' if kwargs.get('cn_plain_fixed_radius') else ''))
            if kwargs.get('cn_plain_fixed_radius'):
                cutoffs = [kwargs.get('cn_plain_fixed_radius')]*len(self.slab_ase)
            else:
                cutoffs = natural_cutoffs(self.slab_ase, mult=1.1)
            nl = NeighborList(cutoffs, skin=0, self_interaction=False, bothways=True)
            nl.update(self.slab_ase)
            cm = nl.get_connectivity_matrix()
            cn_list = [cm[idx].count_nonzero() for idx in surf_sites_indices]


        adsites = []
        adsite_labels = []
        adsites_indices = []
        included_cn_values = []

        max_cn = kwargs.get('max_cn', min(cn_list) + kwargs.get('max_cn_offset'))
        for coords, cn, idx in zip(surf_coords, cn_list, surf_sites_indices):
            if cn <= max_cn:
                adsites.append(coords)
                adsites_indices.append(idx)
                included_cn_values.append(cn)

        if kwargs['selected_sites']: #possibility to select only the main sites, not the surrounding ones
            adsites = [adsites[i] for i in kwargs['selected_sites']]
            adsites_indices = [adsites_indices[i] for i in kwargs['selected_sites']]
            included_cn_values = [included_cn_values[i] for i in kwargs['selected_sites']]

        for coords, cn, idx in zip(adsites, included_cn_values, adsites_indices):
            atom_species = self.asf.slab[idx].species_string
            adsite_labels.append('{0} {1}(cn={2:.1f}),{3:.3f},{4:.3f},'.format(len(adsite_labels), atom_species, cn, *coords[:2]))


        connected_adsites = {}
        all_connected_indices = []
        main_sites_pairs = []
        if kwargs.get('amorphous_surrounding_sites'):
            nn = MinimumDistanceNN(tol=0.2)

            for main_site, label, idx in zip(adsites, adsite_labels, adsites_indices):
                nnsites = nn.get_nn_info(self.asf.slab, idx)
                label_num = int(label.split()[0])
                connected_adsites[label_num] = []
                for nnsite in nnsites:
                    nnindex = nnsite['site_index']
                    if (nnindex in surf_sites_indices if 'surrounding_sites_deltaz' not in kwargs \
                        else abs(main_site[2] - nnsite['site'].coords[2]) < kwargs['surrounding_sites_deltaz']):
                        #print("site: ", f'{label_num}.{len(connected_adsites[label_num])}', "pos:", nnsite['site'].coords)
                        atom_species = self.asf.slab[nnindex].species_string
                        connected_adsites[label_num].append(
                            {
                            'position' : nnsite['site'].coords, 
                            'index' : nnindex,
                            'label' : '{0}.{1} {2},{3:.3f},{4:.3f},'.format(label_num, len(connected_adsites[label_num]), atom_species,  *nnsite['site'].coords[:2]),
                            'duplicate_surrounding' : nnindex in all_connected_indices,
                            'duplicate_main' : {nnindex, idx} in main_sites_pairs #check if already present in swapped order
                            }
                        )
                        #print("duplicate main: ", {nnindex, idx} in main_sites_pairs)
                        all_connected_indices.append(nnindex)
                        if nnindex in adsites_indices:
                            main_sites_pairs.append({idx, nnindex})

        if kwargs.get('VERBOSE'): print('Adsorption sites found.')

        if(kwargs.get('save_image')): #save png to visualize the identified sites
            figname = 'adsorption_sites.png'
            save_adsites_image(
                crystal=False,
                adsites=adsites, 
                adsite_labels=adsite_labels, 
                slab_pymat=self.slab_pymat, 
                connected_adsites = connected_adsites,
                figname=figname, 
                VERBOSE=kwargs.get('VERBOSE'))

        return adsites, adsite_labels, connected_adsites


