import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patheffects as PathEffects
from pymatgen.core import Structure
from pymatgen.analysis.adsorption import plot_slab, get_rot
from ase.visualize.plot import plot_atoms

from xsorb.structures.slab import AdsorptionSite
from xsorb.structures.molecule import MoleculeRotation



def save_rotations_images(mol_rotations_ase : list[MoleculeRotation], 
                          figname : str = 'molecule_orientations.png', 
                          VERBOSE : bool = False):
    '''
    Plot all the rotated molecules from top view

    Args:
    - configs_ase: list of ase Atoms objects containing all the rotated molecules
    - labels: list of strings containing the info on each rotation
    '''

    if VERBOSE: print("Saving image to {0}".format(figname))

    rows_fig = max(int(np.ceil(len(mol_rotations_ase)/5)), 1)
    cols_fig = max(int(np.ceil(len(mol_rotations_ase)/rows_fig)), 1)
    fig = plt.figure(figsize=(3*10 * cols_fig / 5, 3*5 * rows_fig / 3))
    axes = [fig.add_subplot(rows_fig,cols_fig,i) for i in range(1,len(mol_rotations_ase) + 1)]

    for i, rotation in enumerate(mol_rotations_ase):
        conf = rotation.atoms.copy()
        center = (conf.cell[:][0]/2 + conf.cell[:][1]/2)
        conf.translate(center)
        plot_atoms(conf, axes[i], show_unit_cell=2)
        axes[i].set_title('({0}, {1}, {2})'.format(rotation.xrot, rotation.yrot, rotation.zrot))
    fig.suptitle('Molecule orientations (xrot, yrot, zrot)')
    fig.savefig(figname, dpi=800, bbox_inches='tight')

    if VERBOSE: print("Image saved.")



def save_adsites_image(mode : str,
                       adsites : list[AdsorptionSite],        
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

        allowed_modes = ('high_symmetry', 'coord_number', 'coord_number_surrounding')
        if mode not in allowed_modes:
            raise ValueError(f"The mode to plot adsorption sites must be one of {allowed_modes}")
        
        if VERBOSE: print(f"Saving image to {figname}...")
        

        # preparation to plot the sites ####################################################
        fig = plt.figure(figsize=(4,3))
        ax = fig.add_subplot()
        #ax.xaxis.set_tick_params(labelsize=5)
        #ax.yaxis.set_tick_params(labelsize=5)

        #plot slab without the sites, using Pymatgen's function
        plot_slab(slab_pymat, ax, adsorption_sites=False, repeat=3, window=0.7, decay=0.25)


        #calculate useful paraeters to plot the sites
        #w,h = fig.get_size_inches()*fig.dpi
        w = ax.get_xlim()[1] - ax.get_xlim()[0]
        crosses_size = 6.0 * 25. / w
        fontsize     = 2.0 * 25. / w
        mew          = 1.0 * 25. / w
        marker='x'
        sop = get_rot(slab_pymat)   
        ####################################################################################     


        if mode == "high_symmetry":
            ax.set_title('r=ontop, g=bridge, b=hollow')
            cmap_high_sym = {'ontop':'r', 'bridge':'g', 'hollow':'b'}
        else:
            ax.set_title('Adsites based on C.N.')
            from matplotlib import colormaps
            from matplotlib import cm
            from matplotlib.colors import Normalize
            coord_nums = [adsite.coordination_number for adsite in adsites]
            cmap_coord_numb = colormaps.get_cmap('viridis_r')            
            norm = Normalize(vmin=min(coord_nums), vmax=max(coord_nums))
            fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap_coord_numb), ax=ax, label='C.N.') 
         
            
        for main_site in adsites:
            main_site_xy = sop.operate(main_site.coords)[:2]
            
            #in case of coord_number_surrounding, we first draw the surrounding sites
            if mode == "coord_number_surrounding":    
                for surrounding_site in main_site.surrounding_sites:
                    site_xy = sop.operate(surrounding_site.coords)[:2]
                    label = surrounding_site.label
                    if surrounding_site.duplicate_main:
                        label += '^'
                    if surrounding_site.duplicate_surrounding:
                        label += '*'

                    #draw the line connecting the main site to the surrounding site, with point at the end
                    ax.plot([main_site_xy[0], site_xy[0]],
                            [main_site_xy[1], site_xy[1]],
                            '-ok', mfc='r', mec='r', 
                            markersize=crosses_size/3, 
                            linewidth=0.5,
                            mew=mew, 
                            zorder=300000) # zorder to ensure that all crosses are drawn on top
                    #add the label of the surrounding site
                    ax.annotate(label, 
                                xy=site_xy, 
                                xytext=site_xy + np.array([0.2,0 if not surrounding_site.duplicate_surrounding else 0.2]), 
                                fontsize=fontsize*0.8, 
                                path_effects=[PathEffects.withStroke(linewidth=0.3,foreground="w")], 
                                zorder=350000) # zorder to ensure that the text is on top of the crosses

            
            #plot the main sites (or all sites in case of high_symmetry)
            if mode == "high_symmetry":
                color = cmap_high_sym[main_site.type]
            else:
                color = cmap_coord_numb(norm(main_site.coordination_number))
                
            #plot twice to have the edge in black:
            #first plot the cross with a bigger size and black color
            ax.plot(*main_site_xy, 
                    color='black', marker=marker, 
                    markersize=crosses_size*1.05, 
                    mew=mew*1.4, 
                    zorder=500000) # zorder to ensure that all crosses are drawn on top
            #then plot the cross with the right color
            ax.plot(*main_site_xy, 
                    color=color, marker=marker, 
                    markersize=crosses_size, 
                    mew=mew, 
                    zorder=600000) # zorder to ensure that all crosses are drawn on top
            #add the text label on the site
            ax.annotate(str(main_site.label), 
                        xy=main_site_xy, 
                        xytext=main_site_xy, 
                        fontsize=fontsize, 
                        path_effects=[PathEffects.withStroke(linewidth=0.3,foreground="w")], 
                        zorder=1000000) # zorder to ensure that the text is on top of the crosses
        
        fig.savefig(figname, dpi=800, bbox_inches='tight')

        if VERBOSE: print("Image saved.")