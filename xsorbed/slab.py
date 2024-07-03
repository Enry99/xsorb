#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue 28 Feb 2023
@author: Enrico Pedretti

Small helper class to handle the slab

"""

from pymatgen.analysis.adsorption import AdsorbateSiteFinder, plot_slab, get_rot
from pymatgen.analysis.local_env import MinimumDistanceNN
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core import Structure
from ase import Atoms
from ase.io import read
from ase.build.tools import sort
from ase.data import atomic_numbers, covalent_radii
from matplotlib import pyplot as plt
import matplotlib.patheffects as PathEffects
import numpy as np
from ase.constraints import FixCartesian
from ase.geometry import get_distances

from xsorbed import ase_custom


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
        - fixed_indices_slab: list of specific atoms to be fixed 
        (indices start from 0, with the ordering of the atoms in the input file)
        - fix_slab_xyz: which coordinates to fix for the fixed atoms, e.g. [True, True, False] = fix motion in x,y, free to move along z.
        - sort_atoms_by_z: sort the atoms of the slab according to their z coordinate, in reverse order (from higher to lower)
        - translate_slab_from_below_cell_bottom: translate slab to make sure that it is at least 1 Angstrom from the bottom of the cell
    '''

    def __init__(self, slab_filename : str, 
                 layers_threshold : float = 0.5, 
                 surface_sites_height : float = 0.9, 
                 fixed_layers_slab : list = None, 
                 fixed_indices_slab : list = None, 
                 fix_slab_xyz : list = None, 
                 sort_atoms_by_z : bool = False,
                 translate_slab_from_below_cell_bottom : bool = True):
        
        self.slab_ase = ase_custom.Atoms_custom(read(filename=slab_filename))
        self.natoms   = len(self.slab_ase)


        if(min(self.slab_ase.positions[:,2]) < 1 and translate_slab_from_below_cell_bottom):
            #translate slab so that the bottom layer is at least 1 angstrom from the bottom
            self.slab_ase.translate([0,0,1-min(self.slab_ase.positions[:,2])])
        

        #SET CONSTRAINTS (before sorting)##############################
        fixed_atoms_indices = []
        if(fixed_layers_slab): # identify the atoms belonging to the various layers
            fixed_atoms_indices = self._find_atom_indices_fixed_layers(layers_threshold, fixed_layers_slab)
        elif(fixed_indices_slab):
            fixed_atoms_indices = fixed_indices_slab

        c = [FixCartesian(atom_index, mask=[not x for x in fix_slab_xyz]) for atom_index in fixed_atoms_indices]  #we need to negate: in qe 0 = fix, here 1(true)=fix
        self.slab_ase.set_constraint(c) #if no user-defined constraints, c is empty. This is necessary to clean possible constraints read from file (in Molecule there is an if else.)
        ###############################################################

        if sort_atoms_by_z:
            self.slab_ase = sort(self.slab_ase, tags= -self.slab_ase.positions[:, 2])  #sort atoms by height (from higher to lower)

        
        #create pymatgen version of the slab, and initialize the AdsorbateSiteFinder
        
        # First, check if the cell is defined according to the right-handed rule, if not, swap the first two vectors
        # so that the surface normal is found correctly by AdsorbateSiteFinder. Otherwise, the sites will be found on the bottom surface.
        # This also ensures that when using VASP it does not complain.
        cell = self.slab_ase.cell
        if np.dot(np.cross(cell[0], cell[1]), cell[2]) < 0:
            self.slab_ase.set_cell(cell[[1, 0, 2], :])

        import warnings
        with warnings.catch_warnings(): #to suppress the warning about constraints not supported in pymatgen
            warnings.simplefilter("ignore")
            self.slab_pymat = AseAtomsAdaptor.get_structure(self.slab_ase)

        self.asf = AdsorbateSiteFinder(self.slab_pymat, height=surface_sites_height) 


    def _find_atom_indices_fixed_layers(self, layers_threshold : float, fixed_layers_slab : list):
        '''
        Internal helper method to find the indices of the atoms belonging to the 
        layers that we want to fix
        Args:
        - layers_threshold: deltaz for the atoms to be considered as part of the same layer. 
        - fixed_layers_slab: list of layers of the slab to be fixed 
        '''

        #TODO: substitute with ase.geometry.geometry.get_layers
        
        fixed_atoms_indices = []
        original_positions = self.slab_ase.positions
        positions = self.slab_ase.positions.copy().tolist()
        i_layer = 0
        while len(positions): #go on until there are no more atoms
            layer = []
            base_z = min([x[2] for x in positions])

            remove_list = []
            for pos in positions:
                if pos[2] < base_z + layers_threshold:
                    layer.append([idx for idx in range(0, len(original_positions)) if np.allclose(original_positions[idx], pos)][0])
                    remove_list.append(pos)
            for pos in remove_list:
                positions.remove(pos)

            if len(layer) and i_layer in fixed_layers_slab: fixed_atoms_indices += layer
            i_layer += 1

        return fixed_atoms_indices


    def find_adsorption_sites(self, crystal : bool = True, **kwargs):   
        if crystal:
            return self.find_adsorption_sites_crystal(**kwargs)
        else:
            return self.find_adsorption_sites_amorphous(**kwargs)


    def find_adsorption_sites_crystal(self, symm_reduce : float = 0.01, 
                              near_reduce : float = 0.01, 
                              no_obtuse_hollow : bool = True, 
                              selected_sites : list = None,   
                              save_image : bool = False,
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
            sel_adsites = adsites['all'][selected_sites]
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

        if(save_image): #save png to visualize the identified sites
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
            print('Using ase.neighborlist to find coordination numbers.')
            from ase.neighborlist import NeighborList, natural_cutoffs
            cutoffs = natural_cutoffs(self.slab_ase, mult=1.1)
            nl = NeighborList(cutoffs, skin=0, self_interaction=False, bothways=True)
            nl.update(self.slab_ase)
            cm = nl.get_connectivity_matrix()
            cn_list = [cm[idx].count_nonzero() for idx in surf_sites_indices]


        adsites = []
        adsite_labels = []
        adsites_indices = []

        max_cn = kwargs.get('max_cn', min(cn_list) + kwargs.get('max_cn_offset'))
        for coords, cn, idx in zip(surf_coords, cn_list, surf_sites_indices):
            if cn <= max_cn:
                adsites.append(coords)
                atom_species = self.asf.slab[idx].species_string
                adsite_labels.append('{0} {1}(cn={2:.1f}),{3:.3f},{4:.3f},'.format(len(adsite_labels), atom_species, cn, *coords[:2]))
                adsites_indices.append(idx)


        if kwargs['selected_sites']: #possibility to select only the main sites, not the surrounding ones
            adsites = adsites[kwargs['selected_sites']]
            adsite_labels = adsite_labels[kwargs['selected_sites']]
            adsites_indices = adsites_indices[kwargs['selected_sites']]
        

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


    def generate_adsorption_structures(self, molecule : Atoms, 
                                       adsites : list, 
                                       z_distance_from_site : float,
                                       min_z_distance_from_surf : float,
                                       adsites_labels : list,
                                       rotation_label : str,
                                       connected_adsites : dict = None,
                                       surrounding_exclude_main : bool = False,
                                       mol_before_slab : bool = False):
        '''
        Returns the adsorption structures obtained by placing the molecule on all the adsites,
        ensuring that the molecule is not too close to the surface

        Args:
        - molecule: ase Atoms for the molecule to be placed in all the adsorption sites
        - adsites: list of cartesian coordinates of the adsorption sites
        - z_distance_from_site: target distance between reference atom and site.
        - min_z_distance_from_surf: minimum vertical distance between ANY atom of the molecule and ANY atom of the surface
        - adsites_labels: list of labels of the adsorption sites, used to create the final full labels of the adsorption configurations
        - rotation_label: label of the rotation of the molecule, used to create the final full labels of the adsorption configurations
        - mol_before_slab: decide wether to put the molecule before or after the slab in the list of atoms
        '''

        adsorption_structures = []
        full_labels = []

        if not connected_adsites:
            for coords, site_label in zip(adsites, adsites_labels):
                mol = molecule.copy()

                #place the molecule in the adsorption site, then translate it upwards by the target height
                mol.translate(coords) 
                mol.translate([0,0,z_distance_from_site])
                final_deltaz = z_distance_from_site

                #Check for min_z_distance_from_surf and translate accordingly
                dz = mindistance_deltaz(self.slab_ase, mol, min_z_distance_from_surf)
                if dz:
                    mol.translate( [0, 0, dz] )
                    #print("The molecule was translated further by {0} to enforce minimum distance.".format(dz))
                    final_deltaz += dz
                    
                struct = mol + self.slab_ase if mol_before_slab else self.slab_ase + mol
                struct.cell = self.slab_ase.cell
                adsorption_structures.append(struct)
                full_labels.append(rotation_label+site_label+'{:.3f}'.format(final_deltaz))
        
        else:
            yrot = int(float(rotation_label.split(',')[1]))
            for main_site_num, main_site_coords, main_site_label in zip(range(len(adsites)), adsites, adsites_labels):
                
                # create the rotations for when the molecule is horizontal,
                #orienting it towards the nearest neighbors
                
                if yrot != 90:
                    for rel_ads in connected_adsites[main_site_num]:

                        if surrounding_exclude_main and rel_ads['duplicate_main']:
                            continue

                        mol = molecule.copy()

                        #rotate the molecule towards the surrounding site
                        main_to_rel_axis = (rel_ads['position'] - main_site_coords)
                        main_to_rel_axis[2] = 0 #project on xy plane
                        mol.rotate('x', main_to_rel_axis)

                        #place the molecule in the adsorption site, then translate it upwards by the target height
                        mol.translate(main_site_coords) 
                        mol.translate([0,0,z_distance_from_site])
                        final_deltaz = z_distance_from_site

                        #Check for min_z_distance_from_surf and translate accordingly
                        dz = mindistance_deltaz(self.slab_ase, mol, min_z_distance_from_surf)
                        if dz:
                            mol.translate( [0, 0, dz] )
                            #print("The molecule was translated further by {0} to enforce minimum distance.".format(dz))
                            final_deltaz += dz
                            
                        struct = mol + self.slab_ase if mol_before_slab else self.slab_ase + mol
                        struct.cell = self.slab_ase.cell
                        adsorption_structures.append(struct)
                        rotation_label = rotation_label.split(',')
                        rotation_label[2] = 'to_'+rel_ads['label'].split()[0]
                        rotation_label = ','.join(rotation_label)
                        full_labels.append(rotation_label+main_site_label+'{:.3f}'.format(final_deltaz))
                
                else: # vertical molecule
                    rel_ads_positions = [rel_ads['position'] for rel_ads in connected_adsites[main_site_num] \
                                         if not rel_ads['duplicate_surrounding'] and not rel_ads['duplicate_main']]
                    rel_ads_labels = [rel_ads['label'] for rel_ads in connected_adsites[main_site_num] \
                                       if not rel_ads['duplicate_surrounding'] and not rel_ads['duplicate_main']]
                    
                    for coords, site_label in zip([adsites[main_site_num]]+rel_ads_positions, [adsites_labels[main_site_num]]+rel_ads_labels):
                        mol = molecule.copy()

                        #place the molecule in the adsorption site, then translate it upwards by the target height
                        mol.translate(coords) 
                        mol.translate([0,0,z_distance_from_site])
                        final_deltaz = z_distance_from_site

                        #Check for min_z_distance_from_surf and translate accordingly
                        dz = mindistance_deltaz(self.slab_ase, mol, min_z_distance_from_surf)
                        if dz:
                            mol.translate( [0, 0, dz] )
                            #print("The molecule was translated further by {0} to enforce minimum distance.".format(dz))
                            final_deltaz += dz
                            
                        struct = mol + self.slab_ase if mol_before_slab else self.slab_ase + mol
                        struct.cell = self.slab_ase.cell
                        adsorption_structures.append(struct)
                        full_labels.append(rotation_label+site_label+'{:.3f}'.format(final_deltaz))



        return adsorption_structures, full_labels


def closest_pair(slab : Atoms, mol: Atoms):
    '''
    Returns the indices of the closest pair of (slab, molecule) atoms, and their distance

    Args:
    - slab: Atoms object for the slab
    - mol: Atoms object for the molecule 
    '''

    dist_matrix = get_distances(mol.positions, slab.positions, slab.cell, pbc=True)[1]

    i_mol, j_slab = np.unravel_index(np.argmin(dist_matrix), dist_matrix.shape)
    mindist = dist_matrix[i_mol, j_slab]

    return i_mol, j_slab, mindist


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
    while(True):
        #First, find the closest slab-mol atoms pair
        i_mol, j_slab, _ = closest_pair(slab, molcopy)

        #Find the z coordinates of the closest atoms pair, and half their covalent distance
        half_covalent_distance = 0.5 * (covalent_radii[atomic_numbers[slab[j_slab].symbol]] \
                                    + covalent_radii[atomic_numbers[molcopy[i_mol].symbol]])
        zmol = molcopy[i_mol].position[2]
        zslab = slab[j_slab].position[2]

        #Calculate the distance required to enforce the minimum distance
        necessary_min_z_dist = max(min_z_distance_from_surf, half_covalent_distance)
        if(zmol < zslab + necessary_min_z_dist):
            dz = zslab + necessary_min_z_dist - zmol
            dz_tot += dz
            molcopy.translate([0,0,dz])
        else:
            break
    
    return dz_tot


def mol_bonded_to_slab(slab : Atoms, mol: Atoms):
    '''
    Returns true if the molecule is bonded to the slab, i.e. if a
    (slab, molecule) atom pair is closer than 1.1 x sum of covalent radii

    Args:
    - slab: Atoms object for the slab
    - mol: Atoms object for the molecule
    '''

    i_mol, j_slab, distance = closest_pair(slab, mol)

    max_bond_distance = 1.1 * (covalent_radii[atomic_numbers[slab[j_slab].symbol]] \
                                + covalent_radii[atomic_numbers[mol[i_mol].symbol]])
    
    return distance <= max_bond_distance


def save_adsites_image(adsites : list, 
                       adsite_labels : list,
                       slab_pymat : Structure,
                       connected_adsites : dict = None,
                       crystal : bool = True,
                       figname : str = 'adsorption_sites.png',
                       VERBOSE : bool = False):
        '''
        Internal helper function to save an image of the adsorption sites with their numeric label
        Args:
        - adsites: list of cartesian coordinates of the sites
        - adsite_labels: list of the adsites labels (used to color code the sites)
        - slab_pymat: Pymatgen Structure of the slab
        - figname: filename of the image
        '''

        if VERBOSE: print("Saving image to {0}".format(figname))
        

        fig = plt.figure(figsize=(4,3))
        ax = fig.add_subplot()
        #ax.xaxis.set_tick_params(labelsize=5)
        #ax.yaxis.set_tick_params(labelsize=5)

        #plot slab without the sites, using Pymatgen's function
        plot_slab(slab_pymat, ax, adsorption_sites=False, repeat=3, window=0.7, decay=0.25)


        #plot the sites on the slab
        #w,h = fig.get_size_inches()*fig.dpi
        w = ax.get_xlim()[1] - ax.get_xlim()[0]
        crosses_size = 6.0 * 25. / w
        fontsize     = 2.0 * 25. / w
        mew          = 1.0 * 25. / w
        marker='x'

        sop = get_rot(slab_pymat)
        adsites_xy = [sop.operate(ads_site)[:2].tolist() for ads_site in adsites]         

        if not crystal:
            
            if connected_adsites:
                for main_site_idx, related_adsites in connected_adsites.items():     
                    for rel_ads in related_adsites:
                        site_xy = sop.operate(rel_ads['position'])[:2]
                        label = rel_ads['label'].split()[0]
                        if rel_ads['duplicate_main']:
                            label += '^'
                        if rel_ads['duplicate_surrounding']:
                            label += '*'

                        main_site_xy = adsites_xy[main_site_idx]
                        ax.plot([main_site_xy[0], site_xy[0]],
                                [main_site_xy[1], site_xy[1]],
                                '-ok', mfc='r', mec='r', 
                                markersize=crosses_size/3, 
                                linewidth=0.5,
                                mew=mew, 
                                zorder=300000) # zorder to ensure that all crosses are drawn on top
                        ax.annotate(label , 
                                    xy=site_xy, 
                                    xytext=site_xy + np.array([0.2,0 if not rel_ads['duplicate_surrounding'] else 0.2]), 
                                    fontsize=fontsize*0.8, 
                                    path_effects=[PathEffects.withStroke(linewidth=0.3,foreground="w")], 
                                    zorder=350000) # zorder to ensure that the text is on top of the crosses
            
            
            coord_nums = [float(label.split('(')[1].split(')')[0].split('=')[1]) for label in adsite_labels]
            from matplotlib import cm
            from matplotlib.colors import Normalize
            cmap = cm.get_cmap('viridis_r')
            norm = Normalize(vmin=min(coord_nums), vmax=max(coord_nums))

            ax.set_title('Adsites based on C.N.')
            fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, label='C.N.')                        

        else:
            ax.set_title('r=ontop, g=bridge, b=hollow')

       
        for i, label, site_xy in zip(range(len(adsite_labels)), adsite_labels, adsites_xy):
            if crystal:
                if 'ontop' in label:
                    color = 'r'
                elif 'bridge' in label:
                    color = 'g'
                elif 'hollow' in label:
                    color = 'b'
            else:
                color = cmap(norm(coord_nums[i]))
                
            ax.plot(*site_xy, #plot twice to have the edge in black
                    color='black', marker=marker, 
                    markersize=crosses_size*1.05, 
                    mew=mew*1.4, 
                    zorder=500000) # zorder to ensure that all crosses are drawn on top
            ax.plot(*site_xy, 
                    color=color, marker=marker, 
                    markersize=crosses_size, 
                    mew=mew, 
                    zorder=600000) # zorder to ensure that all crosses are drawn on top
            ax.annotate(str(i), 
                        xy=site_xy, 
                        xytext=site_xy, 
                        fontsize=fontsize, 
                        path_effects=[PathEffects.withStroke(linewidth=0.3,foreground="w")], 
                        zorder=1000000) # zorder to ensure that the text is on top of the crosses
        
        fig.savefig(figname, dpi=800, bbox_inches='tight')

        if VERBOSE: print("Image saved.")