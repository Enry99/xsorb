from ase.visualize.plot import plot_atoms
from matplotlib import pyplot as plt
import numpy as np


def save_rotations_images(mol_rotations_ase : list[dict], 
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
        conf = rotation['atoms'].copy()
        center = (conf.cell[:][0]/2 + conf.cell[:][1]/2)
        conf.translate(center)
        plot_atoms(conf, axes[i], show_unit_cell=2)
        axes[i].set_title('({0}, {1}, {2})'.format(* rotation['angles']))
    fig.suptitle('Molecule orientations (xrot, yrot, zrot)')
    fig.savefig(figname, dpi=800, bbox_inches='tight')

    if VERBOSE: print("Image saved.")





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
            from matplotlib import colormaps
            from matplotlib import cm
            from matplotlib.colors import Normalize
            cmap = colormaps.get_cmap('viridis_r')
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