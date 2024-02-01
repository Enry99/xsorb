#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: Enrico Pedretti

Collection of functions for generating images

"""

import numpy as np
import glob, os, shutil
from ase.visualize import view
from ase.io.pov import get_bondpairs
from ase.io import read, write
from xsorbed.slab import Slab
from xsorbed.molecule import Molecule
from xsorbed.settings import Settings
from xsorbed.io_utils import get_calculations_results, _get_configurations_numbers
from xsorbed.dftcode_specific import IN_FILE_PATHS, OUT_FILE_PATHS
from xsorbed.common_definitions import *


#OK (code agnostic)
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


def get_data_for_config_images(calc_type : str, i_or_f = 'f'):
    
    settings = Settings()

    print('Reading files...')

    if i_or_f == 'i': #read from input file
        
        results = None
        calc_indices = _get_configurations_numbers()
        configs = [read(IN_FILE_PATHS[calc_type][settings.program].format(i)) for i in calc_indices]

    elif i_or_f == 'f': #read from output file

        results = get_calculations_results(settings.program, calc_type, settings.E_slab_mol)
        calc_indices = [key for key, val in results['energies'].items() if val is not None]
        configs = [read(OUT_FILE_PATHS[calc_type][settings.program].format(i)) for i in calc_indices]

    
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
    

    #TODO: put here the reindexing for VASP
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


def get_decorated_structures(configs : list, 
                             mol_atoms_indices : list, 
                             slab_atoms_indices : list,
                             ATOM_COLORS_SLAB : list,
                             ATOM_COLORS_MOL : list, 
                             extend_slab : bool = True, 
                             depth_cueing : float = None):
    
    for config in configs:

        mol = config[mol_atoms_indices]
        slab = config[slab_atoms_indices]

        #section to repeat the bulk if mol is partially outside###############
        if(extend_slab):
            slab *= [3,3,1]
            mol.cell = slab.cell          
            mol.translate(+(config.cell[:][0] + config.cell[:][1]) )

            max_a, min_a = ( max(mol.get_scaled_positions()[:,0]), min(mol.get_scaled_positions()[:,0]) )
            max_b, min_b = ( max(mol.get_scaled_positions()[:,1]), min(mol.get_scaled_positions()[:,1]) )
            dx_angstrom = 0.01 #distance in angstrom of surface extending beyond the molecule
            a, b = config.cell.lengths()[:2]
            max_a = max(2/3-0.05/a, max_a + dx_angstrom/a)
            min_a = min(1/3-0.05/a, min_a - dx_angstrom/a)
            max_b = max(2/3-0.05/b, max_b + dx_angstrom/b)
            min_b = min(1/3-0.05/b, min_b - dx_angstrom/b)

            del slab[ [atom.index for atom in slab if (atom.a < min_a or atom.a > max_a or atom.b < min_b or atom.b > max_b)] ]
            Nslab = len(slab)
            slab.translate(-(config.cell[:][0] + config.cell[:][1]) )
            mol.translate(-(config.cell[:][0] + config.cell[:][1]) )
            slab.cell = config.cell
            mol.cell = config.cell
            config = slab + mol
        #######################################################################


    config = configs[0] #just choose one to set the colors
    
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

    return configs, colors, colors_depthcued, textures, textures_depthcued


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
                  extend_slab: bool = True, 
                  depth_cueing : float = None):


    if width_res is None and povray: width_res = 500  # I used 3000. From 1500 is still quite good. 2000 maybe best compromise (still very high res)

    configs, calc_indices, results, ATOM_COLORS_SLAB, ATOM_COLORS_MOL, ATOMIC_RADIUS, BOND_RADIUS, BOND_LINE_WIDTH,\
        CELLLINEWIDTH, mol_atoms_indices, slab_atoms_indices = get_data_for_config_images(calc_type, i_or_f)
    

    configs, colors, colors_depthcued, \
          textures, textures_depthcued = get_decorated_structures(configs=configs,
                                                           mol_atoms_indices=mol_atoms_indices,
                                                           slab_atoms_indices=slab_atoms_indices,
                                                           ATOM_COLORS_SLAB=ATOM_COLORS_SLAB,
                                                           ATOM_COLORS_MOL=ATOM_COLORS_MOL,
                                                           extend_slab=extend_slab,
                                                           depth_cueing=depth_cueing)
    

    if not rotations: #top, lateral
        rot_list = ['', '-5z,-85x'] 
        rotations_labels = ['top', 'lateral']
        color_list = [colors_depthcued, colors]
        textures_list = [textures_depthcued, textures]
    else: #use specified rotation
        rot_list = [rotations]
        rotations_labels = [rotations.replace(',','_')]
        color_list = [colors_depthcued]
        textures_list = [textures_depthcued]


    print('Saving images...')

    main_dir = os.getcwd()
    
    figures_dir = f"{calc_type.lower()}_{images_dirname}/{'initial' if i_or_f == 'i' else 'final'}"    
    os.makedirs(figures_dir, exist_ok=True)
    
    os.chdir(figures_dir)     
    
    for i, config in zip(calc_indices, configs):
        
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
                                        #textures = texture,
                                        camera_dist=1, 
                                        bondatoms=get_bondpairs(config_copy, radius=BOND_RADIUS),
                                        bondlinewidth=BOND_LINE_WIDTH,
                                        #area_light = [(-1.0, -1.0, 200.), 'White', 22.0, 102.0, 20, 2],
                                        )                                
                ).render()
                os.remove(f'{file_label}.pov')
                os.remove(f'{file_label}.ini')

            else:
                write(f'{file_label}.png', config, rotation=rot, scale = 100, colors=color)

    if(i_or_f == 'f'):
        plot_overview_grid(calc_type, rotations_labels[1], calc_indices, povray, results)
        #rotations_labels[0]: use top if not customrot, else use customrot for the grid

    os.chdir(main_dir)

    print(f'All images saved in {figures_dir}.')


#OK (code agnostic)
def view_config(calc_type : str, in_or_out : str, index : int):
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
    elif in_or_out == 'f':
        FILE_PATHS = OUT_FILE_PATHS
    else:
        raise ValueError(f'in or out not recongized. You provided {in_or_out}.')

    file = FILE_PATHS[calc_type][settings.program].format(index)

    try:
        import ase_custom #to make sure that read is correctly monkey-patched
        config = read(file, index=':')
        view(config)
    except:
        print('It was not possible to read the requested configuration.')
            

def relax_animations(povray = False, witdth_res=500, SCREEN_ONLY = False):


    print('Reading files...')
    pwo_list=glob.glob('relax_*.pwo')
    if SCREEN_ONLY:
        configs = [read(file.replace('relax', 'screening'), index=':')[:] for file in pwo_list]
    else:
        configs = [read(file.replace('relax', 'screening'), index=':')[:]+read(file, index=':')[:] for file in pwo_list] #full relax (screening+final)
    labels = [pwo.split('.pwo')[0].split('_')[-1] for pwo in pwo_list]
    print('All files read.')
    
    print('Generating animation(s)...')

    #try:
    settings = Settings()
    slab_filename = settings.slab_filename
    Nbulk = len(read(filename=slab_filename, results_required=False) if slab_filename.split('.')[-1]=='pwo' else read(filename=slab_filename))
 

    from ase.data.colors import jmol_colors
    ATOM_COLORS_SLAB = jmol_colors.copy()
    ATOM_COLORS_MOL  = jmol_colors.copy()

    if os.path.isfile("custom_colors.json"):
        import json
        with open("custom_colors.json", "r") as f:
            custom_colors = json.load(f)

        print("Custom colors read from file.")

        USER_COLORS_SLAB = custom_colors["slab_colors"] if "slab_colors" in custom_colors else []
        USER_COLORS_MOL  = custom_colors["mol_colors"] if "mol_colors" in custom_colors else []
        BOND_RADIUS      = custom_colors["bond_radius"] if "bond_radius" in custom_colors else RADIUS_DEFAULT
        CELLLINEWIDTH     = custom_colors["cell_line_width"] if "cell_line_width" in custom_colors else 0
    else:
        USER_COLORS_SLAB = []
        USER_COLORS_MOL  = []
        BOND_RADIUS      = RADIUS_DEFAULT 
        CELLLINEWIDTH    = 0

    for color in USER_COLORS_SLAB:
        ATOM_COLORS_SLAB[color[0]] = color[1]
    for color in USER_COLORS_MOL:
        ATOM_COLORS_MOL[color[0]] = color[1]

    if(True):
        from ase.build import make_supercell
        Nbulk_original = Nbulk


    if(not os.path.exists('relax_'+images_dirname)):
        os.mkdir('relax_'+images_dirname)
    os.chdir('relax_'+images_dirname)

    if(povray):
        if witdth_res is None: witdth_res = 500 

        for i, config in enumerate(configs): 

            if os.path.isfile('relax_{}_pov.mp4'.format(labels[i])):
                print('relax_{}_pov.mp4 already present. Skipping.'.format(labels[i]))
                continue      
            
            if os.path.exists('temp'):
                shutil.rmtree('temp')
            os.mkdir('temp')  

            os.chdir('temp')

            for j, step in enumerate(config[:]):
                

                #section to repeat the bulk if mol is partially outside###############
                #TODO: prendere le coordinate max tra gli step, in modo da usare sempre quelle per tutti gli step
                if(True):
                    mol = step[Nbulk_original:]
                    slab = step[:Nbulk_original]
                    #print(step.get_scaled_positions())
                    #print(step.cell[:])

                    x_rep, y_rep = (3,3)
                    slab = make_supercell(slab, [[x_rep,0,0], [0,y_rep,0], [0,0,1]], wrap=True) 

                    mol.cell = slab.cell          
                    mol.translate(+(step.cell[:][0] + step.cell[:][1]) )

                    max_a, min_a = ( max(mol.get_scaled_positions()[:,0]), min(mol.get_scaled_positions()[:,0]) )
                    max_b, min_b = ( max(mol.get_scaled_positions()[:,1]), min(mol.get_scaled_positions()[:,1]) )
                    dx_angstrom = 0.1 #distance in angstrom of surface extending beyond the molecule
                    a, b = step.cell.lengths()[:2]
                    max_a = max(2/3-0.05/a, max_a + dx_angstrom/a)
                    min_a = min(1/3-0.05/a, min_a - dx_angstrom/a)
                    max_b = max(2/3-0.05/b, max_b + dx_angstrom/b)
                    min_b = min(1/3-0.05/b, min_b - dx_angstrom/b)

                    del slab[ [atom.index for atom in slab if (atom.a < min_a or atom.a > max_a or atom.b < min_b or atom.b > max_b)] ]
                    Nbulk = len(slab)
                    slab.translate(-(step.cell[:][0] + step.cell[:][1]) )
                    mol.translate(-(step.cell[:][0] + step.cell[:][1]) )
                    slab.cell = step.cell
                    mol.cell = step.cell
                    step = slab + mol
                    #print(step.cell[:])
                #######################################################################
  
                
                colors = [ ATOM_COLORS_SLAB[atom.number] if atom.index < Nbulk else ATOM_COLORS_MOL[atom.number] for atom in step]
                
                step_copy = step.copy()
                step_copy.set_pbc([0,0,0]) #to avoid drawing bonds with invisible replicas
                write('step_{:04d}.pov'.format(j), 
                    step, 
                    format='pov',
                    radii = 0.65, 
                    rotation='-5z,-85x', 
                    povray_settings=dict(canvas_width=witdth_res, celllinewidth=CELLLINEWIDTH, transparent=False, camera_type='orthographic', camera_dist=50., bondatoms=get_bondpairs(step_copy, radius=BOND_RADIUS))
                    #camera_type='perspective'
                ).render()

            os.chdir('..')
            #os.system('convert -delay 20 temp/step_*.png '+'relax_{0}_pov.gif'.format(labels[i]))
            os.system('ffmpeg -framerate 10 -i temp/step_%04d.png  -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2"  -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p '+'relax_{0}_pov.mp4'.format(labels[i]))
            shutil.rmtree('temp')

    else:
        for i, config in enumerate(configs):
            write('relax_{0}.gif'.format(labels[i]), config, rotation='-5z,-85x', interval=150, scale = 100, colors=colors, save_count=None)

    print('All animations saved to {0}.'.format('relax_'+images_dirname))

#OK (code agnostic)
def plot_energy_evolution(calc_type : str):

    from matplotlib import pyplot as plt
    import numpy as np
    
    settings = Settings()

    results = get_calculations_results(program=settings.program, 
                                        calc_type=calc_type,
                                        E_slab_mol=settings.E_slab_mol,
                                        full_evolution=True)

    
    plt.axhline(y=0, linestyle='--', color='black', linewidth=1)
    for i_config, energy_array in results['energies'].items():
        if energy_array[0]:

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