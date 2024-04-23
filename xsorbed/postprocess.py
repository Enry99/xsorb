#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: Enrico Pedretti

Collection of functions for generating images

"""

import numpy as np
import glob, os, shutil
from ase import Atoms
from ase.visualize import view
from ase.io.pov import get_bondpairs
from ase.io import read, write
from xsorbed.slab import Slab
from xsorbed.molecule import Molecule
from xsorbed.settings import Settings
from xsorbed.io_utils import get_calculations_results, _get_configurations_numbers
from xsorbed.dftcode_specific import IN_FILE_PATHS, OUT_FILE_PATHS
from xsorbed.common_definitions import *


def plot_adsorption_sites(ALL : bool = False):
    '''
    Plot an image of the surface with the adsorption sites, with parameters
    according to settings.in

    Args:
    - ALL: if True, plot all sites ignoring symm_reduce  
    '''

    settings = Settings(read_energies=False)

    slab = Slab(slab_filename=settings.slab_filename, 
                surface_sites_height=settings.surface_height, 
                sort_atoms_by_z=settings.sort_atoms_by_z,
                translate_slab_from_below_cell_bottom=settings.translate_slab)

    slab.find_adsorption_sites(symm_reduce = 0 if ALL else settings.symm_reduce, 
                               near_reduce = settings.near_reduce, 
                               selected_sites = None if ALL else settings.selected_sites,
                               save_image = True,
                               figname = 'adsorption_sites_all.png' if ALL else 'adsorption_sites.png',
                               VERBOSE = True)


def get_data_for_config_images(calc_type : str, i_or_f = 'f', read_evolution : bool = False):
    
    settings = Settings()

    print('Reading files...')

    if read_evolution:
        index=':'
    else: index='-1'

    if i_or_f == 'i': #read from input file
        
        results = None
        calc_indices = [i for i in _get_configurations_numbers() if os.path.isfile(IN_FILE_PATHS[calc_type][settings.program].format(i))]
        configs = [read(IN_FILE_PATHS[calc_type][settings.program].format(i)) for i in calc_indices]

    elif i_or_f == 'f': #read from output file

        results = get_calculations_results(settings.program, calc_type, settings.E_slab_mol)
        calc_indices = [key for key, val in results['energies'].items() if val is not None]
        configs = [read(OUT_FILE_PATHS[calc_type][settings.program].format(i), index=index) for i in calc_indices]

    
    #read slab and mol 
    slab = Slab(slab_filename=settings.slab_filename, 
                layers_threshold=settings.layers_height, 
                surface_sites_height=settings.surface_height, 
                fixed_layers_slab=settings.fixed_layers_slab, 
                fixed_indices_slab=settings.fixed_indices_slab, 
                fix_slab_xyz=settings.fix_slab_xyz,
                sort_atoms_by_z=settings.sort_atoms_by_z,
                translate_slab_from_below_cell_bottom=settings.translate_slab)

    mol = Molecule(molecule_filename=settings.molecule_filename,
                   atom_index=settings.selected_atom_index,
                   molecule_axis_atoms=settings.molecule_axis_atoms, 
                   axis_vector=settings.axis_vector, 
                   atoms_subset=settings.mol_subset_atoms, 
                   fixed_indices_mol=settings.fixed_indices_mol, 
                   fix_mol_xyz=settings.fix_mol_xyz)
    

    mol_atoms_indices = np.arange(mol.natoms) if settings.mol_before_slab else np.arange(mol.natoms) + slab.natoms
    slab_atoms_indices = np.arange(slab.natoms) if not settings.mol_before_slab else np.arange(slab.natoms) + mol.natoms

    print('All files read.')

    from ase.data.colors import jmol_colors
    ATOM_COLORS_SLAB = jmol_colors.copy()
    ATOM_COLORS_MOL  = jmol_colors.copy()

    if os.path.isfile(custom_colors_filename):
        import json
        with open(custom_colors_filename, "r") as f:
            custom_colors = json.load(f)

        print("Custom colors read from file.")

        USER_COLORS_SLAB  = custom_colors["slab_colors"] if "slab_colors" in custom_colors else []
        USER_COLORS_MOL   = custom_colors["mol_colors"] if "mol_colors" in custom_colors else []
        ATOMIC_RADIUS     = custom_colors["atomic_radius"] if "atomic_radius" in custom_colors else ATOMIC_RADIUS_DEFAULT
        BOND_RADIUS       = custom_colors["bond_radius"] if "bond_radius" in custom_colors else BOND_RADIUS_DEFAULT
        BOND_LINE_WIDTH   = custom_colors["bond_line_width"] if "bond_line_width" in custom_colors else BOND_LINE_WIDTH_DEFAULT
        CELLLINEWIDTH     = custom_colors["cell_line_width"] if "cell_line_width" in custom_colors else 0
    else:
        USER_COLORS_SLAB  = []
        USER_COLORS_MOL   = []
        ATOMIC_RADIUS     = ATOMIC_RADIUS_DEFAULT
        BOND_RADIUS       = BOND_RADIUS_DEFAULT
        BOND_LINE_WIDTH   = BOND_LINE_WIDTH_DEFAULT
        CELLLINEWIDTH     = 0

    for color in USER_COLORS_SLAB:
        ATOM_COLORS_SLAB[color[0]] = color[1]
    for color in USER_COLORS_MOL:
        ATOM_COLORS_MOL[color[0]] = color[1]


    return configs, calc_indices, results, ATOM_COLORS_SLAB, ATOM_COLORS_MOL, ATOMIC_RADIUS, BOND_RADIUS, BOND_LINE_WIDTH, CELLLINEWIDTH, mol_atoms_indices, slab_atoms_indices


def get_centered_mol_and_slab(mol : Atoms, slab : Atoms):
    """
    Get the centered molecule and slab.

    Parameters:
    - mol (Atoms): The molecule.
    - slab (Atoms): The slab.

    Returns:
    mol (Atoms): The centered molecule.
    slab (Atoms): The centered slab.
    mol_center (ndarray): The center of the molecule in xy, i.e. [xc, yc, 0].
    """

    from ase import neighborlist
    from scipy import sparse

    #avoid modifying the original structure
    mol = mol.copy()
    slab = slab.copy()

    #replicate in both directions to have at least one fully connected molecule
    mol.set_constraint()
    mol_replicated = mol * [2,2,1]
    mol_replicated.set_pbc(False)

    #get the indices of the atoms in each molecule (connected components)
    cutoffs = neighborlist.natural_cutoffs(mol_replicated, mult=1.2)
    neighborList = neighborlist.NeighborList(cutoffs, self_interaction=False, bothways=True)
    neighborList.update(mol_replicated)
    cmatrix = neighborList.get_connectivity_matrix()
    n_components, component_list = sparse.csgraph.connected_components(cmatrix)
    molecules_indices = [ [ atom_id for atom_id in range(len(component_list)) if component_list[atom_id] == molecule_id ] \
                            for molecule_id in range(n_components) ]
    
    #find the molecule with the highest number of atoms (which is therefore fully connected)
    max_molecule_idx = np.argmax([len(molecule) for molecule in molecules_indices])
    mol_replicated = mol_replicated[molecules_indices[max_molecule_idx]]

    #fin the center of mass of the fully connected molecule, that will be translated at the center of the cell
    mol_center = mol_replicated.get_center_of_mass()
    mol_center[2] = 0 #we don't want to translate in the z direction (just xy)

    cell_center = slab.cell[:][0]/2 + slab.cell[:][1]/2  #a1/2 + a2/2 (vector sum)

    mol.translate(cell_center - mol_center)
    slab.translate(cell_center - mol_center)
    #wrap both since we translated them (also the molecule, since we might have chosen a different replica)
    mol.wrap()
    slab.wrap()

    return mol, slab, mol_center


def get_decorated_structures(configs : list, 
                             mol_atoms_indices : list, 
                             slab_atoms_indices : list,
                             ATOM_COLORS_SLAB : list,
                             ATOM_COLORS_MOL : list, 
                             center_molecule : bool = False, 
                             depth_cueing : float = None,
                             cut_vacuum : bool = False):

    cut_height = 3 #cut empty cell to this height above the highest atom
    
    newconfigs = []
    colors_list = []
    colors_depthcued_list = []
    textures_list = []
    textures_depthcued_list = []

    for config in configs:

        if isinstance(config, list): #for full evolution
        
            newconfigs.append([])
            first_frame = True

            for config_j in config:
            
                mol = config_j[mol_atoms_indices]
                slab = config_j[slab_atoms_indices]
                if cut_vacuum: 
                    slab.cell[2][2] = np.max([atom.z for atom in config]) + cut_height
                    if not center_molecule: first_frame = False
                Nslab = len(slab)

                if(center_molecule):
                    #use the same translation for all frame, to avoid jumps
                    if first_frame:
                        _, _, mol_center = get_centered_mol_and_slab(mol, slab)
                        cell_center = config.cell[:][0]/2 + config.cell[:][1]/2  #a1/2 + a2/2
                        first_frame = False

                    mol.translate(cell_center - mol_center)
                    slab.translate(cell_center - mol_center)
                    mol.wrap()
                    slab.wrap()

                newconfigs[-1].append(slab + mol) 
    
        else: #for just single configurations

            mol = config[mol_atoms_indices]
            slab = config[slab_atoms_indices]
            if cut_vacuum: 
                slab.cell[2][2] = np.max([atom.z for atom in config]) + cut_height
            Nslab = len(slab)

            if(center_molecule): 
                mol, slab, _ = get_centered_mol_and_slab(mol, slab)
                
            newconfigs.append(slab + mol)

        
        colors = [ ATOM_COLORS_SLAB[atom.number] if atom.index < Nslab else ATOM_COLORS_MOL[atom.number] for atom in config]
        colors_depthcued = colors.copy()
        textures = ['ase3'] * (len(config))  
        textures_depthcued = textures.copy() 
    

        #fading color for lower layers in top view
        if (depth_cueing):
            zmax = max([atom.z for atom in config if atom.index < Nslab])
            zmin = min([atom.z for atom in config if atom.index < Nslab])
            delta = zmax - zmin
            if depth_cueing < 0:
                raise ValueError("depth_cueing_intensity must be >=0.")
            for atom in config[np.arange(Nslab)]:       
                r,g,b = colors[atom.index] + (np.array([1,1,1]) - colors[atom.index])*(zmax - atom.z)/delta * depth_cueing
                if r>1: r=1
                if g>1: g=1
                if b>1: b=1
                colors_depthcued[atom.index] = [r,g,b]
            
            textures_depthcued = ['pale'] * Nslab + ['ase3'] * (len(config)-Nslab)

        colors_list.append(colors)
        colors_depthcued_list.append(colors_depthcued)
        textures_list.append(textures)
        textures_depthcued_list.append(textures_depthcued) 

    return newconfigs, colors_list, colors_depthcued_list, textures_list, textures_depthcued_list


def plot_overview_grid(calc_type : str, rot_label : str, calc_indices : list, povray : bool, results : dict):

    from matplotlib import pyplot as plt
    import matplotlib.image as mpimg

    #calculate aspect ratio for one image:
    img = mpimg.imread(f"{calc_type.lower()}_{calc_indices[0]}_{rot_label}{'_pov' if povray else ''}.png") 
    ar = float(img.shape[1]) / float(img.shape[0]) #width/height


    n_rows_fig = max(int(np.ceil(len(calc_indices)/5.)), 1)
    n_cols_fig = max(int(np.ceil(len(calc_indices)/float(n_rows_fig))), 1)
    height = n_rows_fig * 2
    width  = n_cols_fig * ar * 2
    fig = plt.figure(figsize=(width, height))
    fig.subplots_adjust(wspace=0.1)
    axes = [fig.add_subplot(n_rows_fig,n_cols_fig,i) for i in range(1,len(calc_indices) + 1)]

    imin = np.argmin([results['energies'][calc_index] for calc_index in calc_indices])
    for i, calc_index in enumerate(calc_indices):
        img = mpimg.imread(f"{calc_type.lower()}_{calc_index}_{rot_label}{'_pov' if povray else ''}.png")
        axes[i].imshow(img)
        axes[i].axis('equal') #ensures all images are in identical rectangular boxes, rescaling the size of the image if necessary
        energy = results['energies'][calc_index]
        status = '**' if results['scf_nonconverged'][calc_index] else '*' if not results['relax_completed'][calc_index] else ''
        axes[i].set_title(f'{energy:.2f}{status} eV', fontsize = 5, pad=1, color='red' if i == imin else 'black')
        #rect = plt.patches.Rectangle((0.0, 0.93), 0.08, 0.07, transform=axes[i].transAxes, facecolor='white', edgecolor='black', linewidth=0.5)
        #axes[i].add_artist(rect)
        #axes[i].text(0.04, 0.985 ...
        axes[i].text(0.045, 0.988, calc_index, fontsize = 3.5, transform=axes[i].transAxes, horizontalalignment="left", verticalalignment="top",
            bbox=dict(boxstyle='square', linewidth=0.5, fc="w", ec="k"),) 
            #fontsize = 4, color="black",horizontalalignment="left", verticalalignment="top")
        axes[i].set_xticks([])
        axes[i].set_yticks([])
    fig.savefig(f"{calc_type.lower()}_overview{'_pov' if povray else ''}.png", dpi=700, bbox_inches='tight')


def config_images(calc_type : str, 
                  i_or_f = 'f', 
                  povray : bool = False, 
                  width_res : int = None, 
                  rotations : str = None, 
                  center_molecule: bool = True, 
                  depth_cueing : float = None,
                  cut_vacuum : bool = False):


    if width_res is None and povray: width_res = 500  # I used 3000. From 1500 is still quite good. 2000 maybe best compromise (still very high res)

    configs, calc_indices, results, ATOM_COLORS_SLAB, ATOM_COLORS_MOL, ATOMIC_RADIUS, BOND_RADIUS, BOND_LINE_WIDTH,\
        CELLLINEWIDTH, mol_atoms_indices, slab_atoms_indices = get_data_for_config_images(calc_type, i_or_f)
    

    if not configs:
        print("First step not completed for any configurations.")
        return
    

    configs, colors_list, colors_depthcued_list, \
          textures_list, textures_depthcued_list = get_decorated_structures(configs=configs,
                                                           mol_atoms_indices=mol_atoms_indices,
                                                           slab_atoms_indices=slab_atoms_indices,
                                                           ATOM_COLORS_SLAB=ATOM_COLORS_SLAB,
                                                           ATOM_COLORS_MOL=ATOM_COLORS_MOL,
                                                           center_molecule=center_molecule,
                                                           depth_cueing=depth_cueing,
                                                           cut_vacuum=cut_vacuum)
    

    print('Saving images...')

    main_dir = os.getcwd()
    
    figures_dir = f"{calc_type.lower()}_{images_dirname}/{'initial' if i_or_f == 'i' else 'final'}"    
    os.makedirs(figures_dir, exist_ok=True)
    
    os.chdir(figures_dir)     

    import xsorbed.ase_custom #monkey patching
    
    for i, config, colors, colors_depthcued, textures, textures_depthcued \
         in zip(calc_indices, configs, colors_list, colors_depthcued_list, textures_list, textures_depthcued_list):


        if not rotations: #lateral, top, 
            rot_list = ['-5z,-85x', ''] 
            rotations_labels = ['lateral', 'top', ]
            color_list = [colors, colors_depthcued]
            textures_list = [textures, textures_depthcued]
        else: #use specified rotation
            rot_list = [rotations]
            rotations_labels = [rotations.replace(',','_')]
            color_list = [colors_depthcued]
            textures_list = [textures_depthcued]
        
        for rot, rot_label, color, texture in zip(rot_list, rotations_labels, color_list, textures_list):

            file_label = f"{calc_type.lower()}_{i}_{rot_label}{'_pov' if povray else ''}"
            
            if(povray): #use POVray renderer (high quality, CPU intensive)
                config_copy = config.copy()
                #config_copy.set_pbc([0,0,0]) #to avoid drawing bonds with invisible replicas

                write(f'{file_label}.pov', 
                    config, 
                    format='pov',
                    radii = ATOMIC_RADIUS, 
                    rotation=rot,
                    colors=color,
                    povray_settings=dict(canvas_width=width_res, 
                                        celllinewidth=CELLLINEWIDTH, 
                                        transparent=False, 
                                        camera_type='orthographic',
                                        textures = texture,
                                        camera_dist=1, 
                                        bondatoms=get_bondpairs(config_copy, radius=BOND_RADIUS),
                                        bondlinewidth=BOND_LINE_WIDTH,
                                        #area_light=[(2., 3., 40.), 'White', .7, .7, 3, 3],
                                        )                                
                ).render()
                os.remove(f'{file_label}.pov')
                os.remove(f'{file_label}.ini')

            else:
                write(f'{file_label}.png', config, rotation=rot, scale = 100, colors=color)

    if(i_or_f == 'f'):
        plot_overview_grid(calc_type, rotations_labels[0], calc_indices, povray, results)
        #rotations_labels[0]: use top if not customrot, else use customrot for the grid

    os.chdir(main_dir)

    print(f'All images saved in {figures_dir}.')



def view_config(calc_type : str, index : int, in_or_out : str):
    '''
    View the selected config with ASE viewer

    Args:
    - calc_type: 'SCREENING' or 'RELAX'
    - in_or_out: read from input or from output
    - index: index of the configuration
    '''

    settings = Settings()
    
    if in_or_out == 'in':
        FILE_PATHS = IN_FILE_PATHS
    elif in_or_out == 'out':
        FILE_PATHS = OUT_FILE_PATHS
    else:
        raise ValueError(f'in or out not recognized. You provided {in_or_out}.')

    file = FILE_PATHS[calc_type][settings.program].format(index)

    try:
        import xsorbed.ase_custom #to make sure that read is correctly monkey-patched
        config = read(file, index=':')
        view(config)
    except:
        print('It was not possible to read the requested configuration.')
            

def relax_animations(calc_type : str,
                     povray : bool = False, 
                     width_res : int = None, 
                     center_molecule: bool = True, 
                     depth_cueing : float = None,
                     cut_vacuum : bool = False):


    if width_res is None and povray: width_res = 500
    

    if calc_type == 'SCREENING':
        configs, calc_indices, results, ATOM_COLORS_SLAB, ATOM_COLORS_MOL, ATOMIC_RADIUS, BOND_RADIUS, BOND_LINE_WIDTH,\
            CELLLINEWIDTH, mol_atoms_indices, slab_atoms_indices = get_data_for_config_images('SCREENING', read_evolution=True)
    else: #screening + relax
        configs_scr, calc_indices_scr, _, _, _, _, _, _, _, _, _ = get_data_for_config_images('SCREENING', read_evolution=True)
        
        configs_rel, calc_indices, results, ATOM_COLORS_SLAB, ATOM_COLORS_MOL, ATOMIC_RADIUS, BOND_RADIUS, BOND_LINE_WIDTH,\
            CELLLINEWIDTH, mol_atoms_indices, slab_atoms_indices = get_data_for_config_images('RELAX', read_evolution=True)

        configs = []

        for conf_rel, i_rel in zip(configs_rel, calc_indices):
            configs.append(configs_scr[calc_indices_scr.index(i_rel)] + conf_rel )

    if not configs:
        print("First step not completed for any configurations.")
        return


    configs, colors_list, colors_depthcued_list, \
        textures_list, textures_depthcued_list = get_decorated_structures(configs=configs,
                                                           mol_atoms_indices=mol_atoms_indices,
                                                           slab_atoms_indices=slab_atoms_indices,
                                                           ATOM_COLORS_SLAB=ATOM_COLORS_SLAB,
                                                           ATOM_COLORS_MOL=ATOM_COLORS_MOL,
                                                           center_molecule=center_molecule,
                                                           depth_cueing=depth_cueing,
                                                           cut_vacuum=cut_vacuum)
    



    print('Generating animation(s)...')

    
    savedir = f'relax_{images_dirname}'
    os.makedirs(savedir, exist_ok=True)

   
    os.chdir(savedir)
    cwd = os.getcwd()

    import xsorbed.ase_custom #monkey patching
    
    for i, config, colors, colors_depthcued, textures, textures_depthcued \
        in zip(calc_indices, configs, colors_list, colors_depthcued_list, textures_list, textures_depthcued_list):


        if(povray):

            if os.path.isfile(f'relax_{config}_pov.mp4'):
                print(f'relax_{config}_pov.mp4 already present. Skipping.')
                continue      
            
            if os.path.exists('temp'):
                shutil.rmtree('temp')
            os.mkdir('temp')  
            os.chdir('temp')

            for j, step in enumerate(config[:]):
                
                step_copy = step.copy()
                #step_copy.set_pbc([0,0,0]) #to avoid drawing bonds with invisible replicas
                write(f'step_{j:04d}.pov', 
                    step_copy, 
                    format='pov',
                    radii = ATOMIC_RADIUS, 
                    colors=colors,
                    rotation='-5z,-85x', 
                    povray_settings=dict(canvas_width=width_res, 
                                        celllinewidth=CELLLINEWIDTH, 
                                        transparent=False, 
                                        camera_type='orthographic',
                                        textures = textures,
                                        camera_dist=1, 
                                        bondatoms=get_bondpairs(step_copy, radius=BOND_RADIUS),
                                        bondlinewidth=BOND_LINE_WIDTH,
                                        #area_light=[(2., 3., 40.), 'White', .7, .7, 3, 3],
                                        )   
                ).render()

            os.chdir(cwd)
            #os.system('convert -delay 20 temp/step_*.png '+'relax_{0}_pov.gif'.format(labels[i]))
            os.system(f'ffmpeg -framerate 8 -i temp/step_%04d.png  -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" \
                       -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p relax_{i}_pov.mp4')
            shutil.rmtree('temp')

        else:
            write(f'relax_{i}.gif', config, rotation='-5z,-85x', interval=150, scale = 100, colors=colors, save_count=None)


    print(f'All animations saved to {savedir}.')


def plot_energy_evolution(calc_type : str):

    from matplotlib import pyplot as plt
    import numpy as np
    
    settings = Settings()

    results = get_calculations_results(program=settings.program, 
                                        calc_type=calc_type,
                                        E_slab_mol=settings.E_slab_mol,
                                        full_evolution=True)

    
    if 0 not in settings.E_slab_mol:
        plt.axhline(y=0, linestyle='--', color='black', linewidth=1)
    for i_config, energy_array in results['energies'].items():
        if energy_array is not None and len(energy_array) > 0:

            #completion status of the calculations
            if results['relax_completed'][i_config]:
                status = ''
            elif results['scf_nonconverged'][i_config]:
                status = '**'
            else: 
                status = '*' #relax not completed, but scf converged
            
            plt.plot(energy_array, '-', label=f'{i_config}: {energy_array[-1]:.2f}{status} eV')

            if ('*' in status): 
                symbol = 'x' if status == '*' else '^'
                color = 'black' if status == '*' else 'red'
                plt.plot(len(energy_array)-1, energy_array[-1], symbol, color=color)
        else:
            print(f'Config. {i_config} job has not reached the first scf convergence. It will be skipped.')
            
    
    plt.title('Energy evolution during optimization')
    plt.xlabel('step')
    plt.ylabel('energy (eV)')
    plt.grid(linestyle='dotted')
    plt.legend(title="Config, energy", 
               ncols=np.ceil(len([x for x in results['energies'].values() if x is not None])/10), 
               prop={'size': 6  if calc_type == 'SCREENING' else 8})
    energy_plot_filename = f'{calc_type.lower()}_energies.png'
    plt.savefig(energy_plot_filename, dpi=300, bbox_inches='tight')
    print(f'Plot saved in {energy_plot_filename}')