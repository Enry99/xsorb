#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue 28 Feb 2023
@author: Enrico Pedretti

Small helper class to deal with the slab

"""

from pymatgen.analysis.adsorption import AdsorbateSiteFinder, plot_slab, get_rot
from pymatgen.analysis.local_env import MinimumDistanceNN
from pymatgen.io import ase
from ase.io import read
from ase.build.tools import sort
from matplotlib import pyplot as plt
import numpy as np


class Slab:

    def __init__(self, slab_filename : str, threshold = 0.2):
        '''
            Read slab from file (e.g. Quantum ESPRESSO pwi/pwo or .xyz)
        '''
        self.slab_ase = read(filename=slab_filename, results_required=False) if slab_filename.split('.')[-1]=='pwo' else read(filename=slab_filename)
        self.slab_ase.set_initial_magnetic_moments(self.slab_ase.get_global_number_of_atoms()*[0])

        #translate slab so that the bottom layer is at least 1 angstrom from the bottom
        zmin = min(self.slab_ase.positions[:,2])
        if(zmin < 1):
            self.slab_ase.translate([0,0,1-zmin])
        

        #identification of the layers##################################
        slab = sort(self.slab_ase, tags= self.slab_ase.positions[:, 2])
        self.layers = [[]] #each list element corresponds to a layer, and it is itself a list of the
        #atoms contained in the layer
        zmin = min(slab.positions[:,2])
        i_layer = 0
        for i, z in enumerate(slab.positions[:,2]):
            if(z-zmin > threshold): #new layer
                zmin = z
                i_layer = i_layer+1
                self.layers.append([])
            self.layers[i_layer].append(i)
        ###############################################################

        #reindex mapping before sorting by z for fixing atoms by index#
        self.reindex_map = np.argsort(-self.slab_ase.positions[:, 2], kind='stable')
        #NOTE: here we used np.argsort to get the indices, while ase.sort (which uses sorted library) to actually sort the elements.
        # Check if same result (possible problems with very similar numbers)         
        ###############################################################      
    
        self.slab_ase = sort(self.slab_ase, tags= -self.slab_ase.positions[:, 2])  #sort atoms by height (from higher to lower)
        self.slab_pymat = ase.AseAtomsAdaptor.get_structure(self.slab_ase)
        self.asf = AdsorbateSiteFinder(self.slab_pymat)

    def get_atoms_by_layers(self, layers_list : list[int]):
        atoms_list = []
        for layer in layers_list:
            atoms_list += self.layers[layer]
        return atoms_list


    def find_adsorption_sites(self, distance_from_surf=0., symm_reduce_thr=0.01, near_reduce_thr=0.01, no_obtuse_hollow=True, save_image = False, selected_sites : list = []):
        '''
        Returns a list of cartesian coordinates of the adsites, and a list of labels ('ontop', x, y).
        Optionally it saves a figure with the sites on the surface.

        Args:
            distance_from_surf: distance of the site (where the selected atom of the molecule is to be placed) from the surface
            symm_reduce_thr, near_reduce_thr: thresholds for removing sites duplicates (increase them to reduce duplicates).
            no_obtuse_hollow: avoid considering hollow sites inside obtuse triangles of the Delaunay triangulation of topmost layer used to find sites.
        '''

        if symm_reduce_thr == 0:
            figname = 'adsorption_sites_all.png'
        else:
            figname = 'adsorption_sites.png'

        adsites = self.asf.find_adsorption_sites(distance=distance_from_surf, symm_reduce=symm_reduce_thr, near_reduce=near_reduce_thr, no_obtuse_hollow=no_obtuse_hollow)
        if selected_sites:
            sel_adsites = [adsites['all'][i] for i in selected_sites]
        else:
            sel_adsites = adsites['all']
        

        adsite_labels = []
        nn = MinimumDistanceNN()
        #run over all the slab_adsites, classifiying them by checking if the
        #i-th element of 'all' is in one of the three lists 'ontop', 'hollow', 'bridge'
        for site in sel_adsites:
            #dummy structure just to place one atom in the site
            slab = self.slab_pymat.copy()
            coords = site.tolist()
            coords[2] += 2
            slab.append('O', coords, coords_are_cartesian=True)
            coord_n = nn.get_cn(slab, len(slab.sites)-1)
            nn_list = nn.get_nn(slab, len(slab.sites)-1)       

            if any((site == x).all() for x in adsites['ontop']):
                first_nn_species = nn_list[0].species_string
                adsite_labels.append('ontop_{0},{1:.3f},{2:.3f},'.format(first_nn_species, *site[:2]))   
            elif any((site == x).all() for x in adsites['hollow']):
                adsite_labels.append('hollow_c{0},{1:.3f},{2:.3f},'.format(coord_n, *site[:2])) 
            else:
                if(coord_n>=4): #attemps to fix the problem of fake bridges for 4-fold sites
                    adsite_labels.append('hollow_c{0},{1:.3f},{2:.3f},'.format(coord_n, *site[:2]))
                else:
                    if len(nn_list) >=2:
                        distance = np.linalg.norm(nn_list[0].coords[:2] - nn_list[1].coords[:2])
                        adsite_labels.append('bridge_{0:.2f},{1:.3f},{2:.3f},'.format(distance, *site[:2]))
                    else: adsite_labels.append('bridge,{0:.3f},{1:.3f},'.format(*site[:2]))
        
        if(save_image): #save png to visualize the identified sites
            sop = get_rot(self.asf.slab)

            fig = plt.figure(figsize=(4,3))
            ax = fig.add_subplot(111)
            plot_slab(self.slab_pymat, ax, adsorption_sites=False, window=0.7)

            adsites_xy = [sop.operate(ads_site)[:2].tolist() for ads_site in sel_adsites]
            for i, site in enumerate(sel_adsites):
                if any((site == x).all() for x in adsites['ontop']):
                    color = 'r'
                elif any((site == x).all() for x in adsites['bridge']):
                    color = 'g'
                elif any((site == x).all() for x in adsites['hollow']):
                    color = 'b'
                ax.plot(*adsites_xy[i], color=color, marker="x", markersize=3, mew=1, linestyle="", zorder=500000)
                ax.annotate(str(i), xy=adsites_xy[i], xytext=adsites_xy[i], fontsize=2, zorder=1000000)
                            
            ax.set_title('Adsites: r=ontop, g=bridge, b=hollow')
            fig.savefig(figname, dpi=1500, bbox_inches='tight')

        return sel_adsites, adsite_labels

    def generate_adsorption_structures(self, molecule, adsites):
        
        structs = []

        for coords in adsites:
            mol = molecule.copy()
            mol.translate(coords)
            structs.append(self.slab_ase + mol)

        return structs


#NOTE!!!: CURRENTLY NOT WORKING
def adsorb_both_surfaces(all_mol_on_slab_configs_pymat : list): 
    new_adslabs = []
    
    for adslab in all_mol_on_slab_configs_pymat: 

        # Find the adsorbate sites and indices in each slab           
        _, adsorbates, indices = False, [], []
        for i, site in enumerate(adslab.sites):
            if site.surface_properties == "adsorbate":
                    adsorbates.append(site)
                    indices.append(i)

        # Clean slab to find center of mass
        slab_noads = adslab.copy()
        slab_noads.remove_sites(indices)

        new_slab = adslab.copy()
        for adsorbate in adsorbates:
            p2 = 0 #qui ci vogliono le coordinate di asdorbate con la z invertita
            #rispetto al centro geometrico (non di massa) della slab_noads
            new_slab.append(adsorbate.specie, p2, properties={"surface_properties": "adsorbate"})
        new_adslabs.append(new_slab)
    
    return new_adslabs

