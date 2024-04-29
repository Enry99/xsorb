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
from ase.build import make_supercell

from xsorbed import ase_custom

#NOTE: ALL MODIFICATIONS DONE
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
        slabcopy = self.slab_ase.copy() 
        del slabcopy.constraints #to suppress the warning about constraints not supported in pymatgen
        self.slab_pymat = AseAtomsAdaptor.get_structure(slabcopy)
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


    def find_adsorption_sites(self, symm_reduce : float = 0.01, 
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

        if(save_image): #save png to visualize the identified sites
            save_adsites_image(sel_adsites, adsite_labels, self.slab_pymat, figname, VERBOSE)

        return sel_adsites, adsite_labels


    def generate_adsorption_structures(self, molecule : Atoms, 
                                       adsites : list, 
                                       z_distance_from_site : float,
                                       min_z_distance_from_surf : float,
                                       adsites_labels : list,
                                       rotation_label : str,
                                       mol_before_slab = False):
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
                

            adsorption_structures.append(mol + self.slab_ase if mol_before_slab else self.slab_ase + mol)
            full_labels.append(rotation_label+site_label+'{:.3f}'.format(final_deltaz))

        return adsorption_structures, full_labels


def closest_pair(slab : Atoms, mol: Atoms):
    '''
    Returns the indices of the closest pair of (slab, molecule) atoms, and their distance

    Args:
    - slab: Atoms object for the slab
    - mol: Atoms object for the molecule 
    '''
    distances = []
    for mol_atom in mol:
        for slab_atom in slab:
            distances.append(np.linalg.norm(mol_atom.position - slab_atom.position))
    mindist = min(distances)

    i_min = distances.index(mindist)
    i_mol, j_slab = ( i_min // len(slab), i_min % len(slab) )

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

    #Make replicas to conisder the pbcs
    #TODO: in the future, change this with neighbours, instead of replicating the cell
    slab_rep = make_supercell(slab, [[3,0,0], [0,3,0], [0,0,1]], wrap=False)
    point = slab.cell[:][0] + slab.cell[:][1]
    slab_rep.translate(-point) 

    dz_tot = 0
    molcopy = mol.copy()
    while(True):
        #First, find the closest slab-mol atoms pair
        i_mol, j_slab, _ = closest_pair(slab_rep, molcopy)

        #Find the z coordinates of the closest atoms pair, and half their covalent distance
        half_covalent_distance = 0.5 * (covalent_radii[atomic_numbers[slab_rep[j_slab].symbol]] \
                                    + covalent_radii[atomic_numbers[molcopy[i_mol].symbol]])
        zmol = molcopy[i_mol].position[2]
        zslab = slab_rep[j_slab].position[2]

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
        ax = fig.add_subplot(111)
        #ax.xaxis.set_tick_params(labelsize=5)
        #ax.yaxis.set_tick_params(labelsize=5)

        #plot slab without the sites, using Pymatgen's function
        plot_slab(slab_pymat, ax, adsorption_sites=False, window=1.0, decay=0.25)


        #plot the sites on the slab
        #w,h = fig.get_size_inches()*fig.dpi
        w = ax.get_xlim()[1] - ax.get_xlim()[0]
        crosses_size = 6.0 * 25. / w
        fontsize     = 2.0 * 25. / w
        mew          = 1.0 * 25. / w

        sop = get_rot(slab_pymat)
        adsites_xy = [sop.operate(ads_site)[:2].tolist() for ads_site in adsites]
        for i, label, site_xy in zip(range(len(adsite_labels)), adsite_labels, adsites_xy):
            if 'ontop' in label:
                color = 'r'
            elif 'bridge' in label:
                color = 'g'
            elif 'hollow' in label:
                color = 'b'
            ax.plot(*site_xy, 
                    color=color, marker="x", 
                    markersize=crosses_size, 
                    mew=mew, 
                    linestyle="", 
                    zorder=500000) # zorder to ensure that all crosses are drawn on top
            ax.annotate(str(i), 
                        xy=site_xy, 
                        xytext=site_xy, 
                        fontsize=fontsize, 
                        path_effects=[PathEffects.withStroke(linewidth=0.25,foreground="w")], 
                        zorder=1000000) # zorder to ensure that the text is on top of the crosses
                        
        ax.set_title('r=ontop, g=bridge, b=hollow')
        fig.savefig(figname, dpi=800, bbox_inches='tight')

        if VERBOSE: print("Image saved.")