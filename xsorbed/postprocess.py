#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: Enrico Pedretti

Collection of functions for generating images

"""

from ase.visualize import view
from ase.io.pov import get_bondpairs
from ase.io import read, write
import glob, sys, os, shutil
from natsort import natsorted
from slab import Slab
from settings import Settings
from filenames import *




def read_energy(filename: str, Eslab = 0, Emol = 0):
    with open(filename, 'r') as f:
        pwo = f.readlines()

    toten = 0
    scf_terminated = False
    for line in pwo: #make sure to get the last one (useful in relaxations)
        if '!' in line: 
            toten = line.split()[4]
        if 'End of self-consistent calculation' in line:
            scf_terminated = True
    if(scf_terminated and toten != 0): 
        return (float(toten) - Eslab - Emol)*rydbergtoev
    else: return None


def plot_adsorption_sites(ALL = False):

    settings = Settings(read_energies=False)
    slab = Slab(settings.slab_filename, surface_sites_height=settings.surface_height)
  
    slab.find_adsorption_sites(* {
        "distance":0, 
        'symm_reduce': 0 if ALL else settings.symm_reduce, 
        'near_reduce':settings.near_reduce, 
        'no_obtuse_hollow':True}.values(),
        selected_sites= [] if ALL else settings.selected_sites,
        ALL = ALL,
        save_image=True
        )


def config_images(which : str, i_or_f = 'f', povray = False, witdth_res=500, index : str = None, rotations : str = None):

    if witdth_res is None and povray: witdth_res = 500  # I used 3000. From 1500 is still quite good. 2000 maybe best compromise (still very high res)
    if i_or_f == 'i':
        pw = 'pwi'
    else:
        pw = 'pwo'
    print('Reading files...')
    pw_list=natsorted(glob.glob(pw_files_prefix+which+'_*.'+pw))
    if(not pw_list):
        print("Files not found. Quitting.")
        sys.exit(1)

    configs = []
    labels = []
    uncompleted = []
    for file in pw_list:
        with open(file, 'r') as f:
            lines = f.readlines()
        STOP = False
        for line in lines:
            if 'convergence NOT achieved' in line:
                STOP = True            
            if 'bfgs converged' in line:
                STOP = False
                break          
        if(STOP): 
            print("Found 'convergence NOT achieved' in {0}, so relaxation was not completed. It will be skipped.".format(file))
            uncompleted.append(file)
            continue
        configs.append(read(file) if pw == 'pwi' else read(file, results_required=False))
        labels.append(file.split('.'+pw)[0].split('_')[-1])

    print('All files read.')

    if index is not None:
        configs = [configs[labels.index(index)]] #just one
        labels = [labels[labels.index(index)]]

    print('Saving images...')

    #try:
    settings = Settings()
    E_slab_mol = settings.E_slab_mol
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
    else:
        USER_COLORS_SLAB = []
        USER_COLORS_MOL  = []
        BOND_RADIUS      = RADIUS_DEFAULT

    for color in USER_COLORS_SLAB:
        ATOM_COLORS_SLAB[color[0]] = color[1]
    for color in USER_COLORS_MOL:
        ATOM_COLORS_MOL[color[0]] = color[1]

    if(True):
        from ase.build import make_supercell
        Nbulk_original = Nbulk

    main_dir = os.getcwd()

    figures_dir = '{0}/{1}'.format(which+'_'+images_dirname, 'initial' if i_or_f == 'i' else 'final')

    if(not os.path.exists(figures_dir)):
        os.makedirs(figures_dir)
    os.chdir(figures_dir) 

    for i, config in enumerate(configs):
        label = labels[i]

        #section to repeat the bulk if mol is partially outside###############
        if(True):
            mol = config[Nbulk_original:]
            slab = config[:Nbulk_original]
            #print(config.get_scaled_positions())
            #print(config.cell[:])

            x_rep, y_rep = (3,3)
            slab = make_supercell(slab, [[x_rep,0,0], [0,y_rep,0], [0,0,1]], wrap=True) 

            mol.cell = slab.cell          
            mol.translate(+(config.cell[:][0] + config.cell[:][1]) )

            max_a, min_a = ( max(mol.get_scaled_positions()[:,0]), min(mol.get_scaled_positions()[:,0]) )
            max_b, min_b = ( max(mol.get_scaled_positions()[:,1]), min(mol.get_scaled_positions()[:,1]) )
            dx_angstrom = 0.1 #distance in angstrom of surface extending beyond the molecule
            a, b = config.cell.lengths()[:2]
            max_a = max(2/3-0.01/a, max_a + dx_angstrom/a)
            min_a = min(1/3-0.01/a, min_a - dx_angstrom/a)
            max_b = max(2/3-0.01/b, max_b + dx_angstrom/b)
            min_b = min(1/3-0.01/b, min_b - dx_angstrom/b)

            del slab[ [atom.index for atom in slab if (atom.a < min_a or atom.a > max_a or atom.b < min_b or atom.b > max_b)] ]
            Nbulk = len(slab)
            slab.translate(-(config.cell[:][0] + config.cell[:][1]) )
            mol.translate(-(config.cell[:][0] + config.cell[:][1]) )
            slab.cell = config.cell
            mol.cell = config.cell
            config = slab + mol
            #print(config.cell[:])
        #######################################################################


        #fading color for lower layers in top view
        zmax = max([atom.z for atom in config if atom.index < Nbulk])
        zmin = min([atom.z for atom in config if atom.index < Nbulk])
        delta = zmax - zmin
        import numpy as np

        colors = [ ATOM_COLORS_SLAB[atom.number] if atom.index < Nbulk else ATOM_COLORS_MOL[atom.number] for atom in config]
        colors_top = [ (colors[atom.index] + (np.array([1,1,1]) - colors[atom.index])*(zmax - atom.z)/delta).round(4) if atom.index < Nbulk else colors[atom.index] for atom in config ]

        if(povray):
            config_copy = config.copy()
            config_copy.set_pbc([0,0,0]) #to avoid drawing bonds with invisible replicas

            if rotations is not None: #use specified rotation
                write(pw_files_prefix+which+'_{0}_{1}_pov.pov'.format(label, rotations.replace(',','_')), 
                    config, 
                    format='pov',
                    radii = 0.65, 
                    rotation=rotations,
                    colors=colors,
                    povray_settings=dict(canvas_width=witdth_res, celllinewidth=0, transparent=False, camera_type='orthographic', camera_dist=50., bondatoms=get_bondpairs(config_copy, radius=BOND_RADIUS))
                    #camera_type='perspective'
                ).render()
                os.remove(pw_files_prefix+which+'_{0}_{1}_pov.pov'.format(label, rotations.replace(',','_')))
                os.remove(pw_files_prefix+which+'_{0}_{1}_pov.ini'.format(label, rotations.replace(',','_')))
            else: #front and top view
                #front view
                write(pw_files_prefix+which+'_{0}_pov.pov'.format(label), 
                    config, 
                    format='pov',
                    radii = 0.65, 
                    rotation='-5z,-85x', 
                    colors=colors,
                    povray_settings=dict(canvas_width=witdth_res, celllinewidth=0, transparent=False, camera_type='orthographic', camera_dist=50., bondatoms=get_bondpairs(config_copy, radius=BOND_RADIUS))
                    #camera_type='perspective'
                ).render()
                os.remove(pw_files_prefix+which+'_{0}_pov.pov'.format(label))
                os.remove(pw_files_prefix+which+'_{0}_pov.ini'.format(label))

                #top view

                textures = ['pale'] * Nbulk + ['ase3'] * (len(config)-Nbulk)
                write(pw_files_prefix+which+'_{0}_top_pov.pov'.format(label), 
                    config, 
                    format='pov',
                    radii = 0.65, 
                    colors=colors_top,
                    povray_settings=dict(canvas_width=witdth_res, celllinewidth=0, transparent=False, textures = textures,
                        camera_type='orthographic', camera_dist=50., bondatoms=get_bondpairs(config_copy, radius=BOND_RADIUS))
                    #camera_type='perspective'
                ).render()
                os.remove(pw_files_prefix+which+'_{0}_top_pov.pov'.format(label))
                os.remove(pw_files_prefix+which+'_{0}_top_pov.ini'.format(label))

        else:
            if rotations is not None: #use specified rotation
                write(pw_files_prefix+which+'_{0}_{1}.png'.format(label, rotations), config, rotation=rotations, scale = 100, colors=colors)
            else: #front and top view           
                write(pw_files_prefix+which+'_{0}.png'.format(label), config, rotation='-10z,-80x', scale = 100, colors=colors) #front view
                write(pw_files_prefix+which+'_{0}_top.png'.format(label), config, scale = 100, colors=colors_top) #top view

    os.chdir(main_dir)

    if(pw=='pwo' and index is None):
        
        import numpy as np
        from matplotlib import pyplot as plt
        import matplotlib.image as mpimg

        #calculate aspect ratio for one image:
        img = mpimg.imread(figures_dir+'/'+pw_files_prefix+which+'_{0}{1}.png'.format(labels[0], '_pov' if povray else ''))
        ar = float(img.shape[1]) / float(img.shape[0]) #width/height

        n_rows_fig = max(int(np.ceil(len(configs)/5.)), 1)
        n_cols_fig = max(int(np.ceil(len(configs)/float(n_rows_fig))), 1)
        height = n_rows_fig * 2
        width  = n_cols_fig * ar * 2
        fig = plt.figure(figsize=(width, height))
        fig.subplots_adjust(wspace=0.1)
        axes = [fig.add_subplot(n_rows_fig,n_cols_fig,i) for i in range(1,len(configs) + 1)]
    
        energies = [read_energy(file, *E_slab_mol) for file in pw_list if file not in uncompleted]

        for i, conf in enumerate(configs):
            label = labels[i]
            img = mpimg.imread(figures_dir+'/'+pw_files_prefix+which+'_{0}{1}.png'.format(label, '_pov' if povray else ''))
            axes[i].imshow(img)
            axes[i].axis('equal') #ensures all images are in identical rectangular boxes, rescaling the size of the image if necessary
            axes[i].set_title('{0:.2f} eV'.format(energies[i]), fontsize = 5, pad=1)
            #rect = plt.patches.Rectangle((0.0, 0.93), 0.08, 0.07, transform=axes[i].transAxes, facecolor='white', edgecolor='black', linewidth=0.5)
            #axes[i].add_artist(rect)
            #axes[i].text(0.04, 0.985 ...
            axes[i].text(0.045, 0.988, label, fontsize = 3.5, transform=axes[i].transAxes, horizontalalignment="left", verticalalignment="top",
                bbox=dict(boxstyle='square', linewidth=0.5, fc="w", ec="k"),) 
                #fontsize = 4, color="black",horizontalalignment="left", verticalalignment="top")
            axes[i].set_xticks([])
            axes[i].set_yticks([])
        fig.savefig('{1}/{0}_overview{2}.png'.format(which, figures_dir,'_pov' if povray else ''), dpi=700, bbox_inches='tight')

    print('All images saved in {0}.'.format(which+'_'+images_dirname))


def view_config(which : str, index : int):
    if which == 's':
        which = 'screening'
        pw = 'pwi'
    elif which == 'r':
        which = 'relax'
        pw = 'pwo'
    file = pw_files_prefix+which+'_{0}.'.format(index)+pw
    try:
        config = read(file, index=':')
    except:
        config =  read(file) if pw == 'pwi' else read(file, results_required=False)
        print('Not possible to load the full relaxation history because of the ase bug for Espresso > 6.7. Showing only the last configuration.')
    view(config)        


def relax_animations(povray = False, witdth_res=500, SCREEN_ONLY = False):


    print('Reading files...')
    pwo_list=glob.glob(pw_files_prefix+'relax_*.pwo')
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
    else:
        USER_COLORS_SLAB = []
        USER_COLORS_MOL  = []
        BOND_RADIUS      = RADIUS_DEFAULT

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
                    max_a = max(2/3-0.01/a, max_a + dx_angstrom/a)
                    min_a = min(1/3-0.01/a, min_a - dx_angstrom/a)
                    max_b = max(2/3-0.01/b, max_b + dx_angstrom/b)
                    min_b = min(1/3-0.01/b, min_b - dx_angstrom/b)

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
                write(pw_files_prefix+'step_{:04d}.pov'.format(j), 
                    step, 
                    format='pov',
                    radii = 0.65, 
                    rotation='-5z,-85x', 
                    povray_settings=dict(canvas_width=witdth_res, celllinewidth=0, transparent=False, camera_type='orthographic', camera_dist=50., bondatoms=get_bondpairs(step_copy, radius=BOND_RADIUS))
                    #camera_type='perspective'
                ).render()

            os.chdir('..')
            #os.system('convert -delay 20 temp/step_*.png '+pw_files_prefix+'relax_{0}_pov.gif'.format(labels[i]))
            os.system('ffmpeg -framerate 10 -i temp/step_%04d.png  -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2"  -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p '+pw_files_prefix+'relax_{0}_pov.mp4'.format(labels[i]))
            shutil.rmtree('temp')

    else:
        for i, config in enumerate(configs):
            write(pw_files_prefix+'relax_{0}.gif'.format(labels[i]), config, rotation='-5z,-85x', interval=150, scale = 100, colors=colors, save_count=None)

    print('All animations saved to {0}.'.format('relax_'+images_dirname))


def plot_energy_evolution(which='relax'):

    settings = Settings()
    Eslab, Emol = (settings.E_slab_mol[0], settings.E_slab_mol[1])

    print('Reading files...')
    pwo_list=natsorted(glob.glob(pw_files_prefix+which+'_*.pwo'))
    labels = [int(pwo.split('.pwo')[0].split('_')[-1]) for pwo in pwo_list]

    totens = []
    relax_terminated = []
    for file in pwo_list:

        end = False
        totens.append([])

        with open(file, 'r') as f:
            pwo = f.readlines()

        for line in pwo: #make sure to get the last one (useful in relaxations)
            if '!' in line: 
                totens[-1].append( (float(line.split()[4]) - (Eslab+Emol)) * rydbergtoev )
            if 'Final energy' in line:
                end = True
        relax_terminated.append(end)
    print('All files read.')


    from matplotlib import pyplot as plt
    plt.axhline(y=0, linestyle='--', color='black', linewidth=1)
    for i, config_e in enumerate(totens):
        if(False):#len(config_e) > 10): #skip first 10 steps
            plt.plot([*range(10, len(config_e))], config_e[10:], '-', label=labels[i])
            plt.xlim(xmin=10)
        else:
            if config_e:
                plt.plot(config_e, '-', label='{0}: {1:.2f}{2} eV'.format(labels[i], config_e[-1], '' if relax_terminated[i] else '*'))
                if (not relax_terminated[i]): plt.plot(len(config_e)-1, config_e[-1], 'x', color='black')
            else:
                print('Config. {0} job has not reached scf convergence. It will be skipped.'.format(labels[i]))
            #plt.xlim(xmin=0)
            
    
    plt.title('Relax energies')
    plt.xlabel('step')
    plt.ylabel('energy (eV)')
    plt.grid(linestyle='dotted')
    import math
    plt.legend(title="Config, energy", ncols=math.ceil(len(totens)/10), prop={'size': 6  if which == 'screening' else 8})
    energy_plot_filename = '{0}_energies.png'.format(which)
    plt.savefig(energy_plot_filename, dpi=300, bbox_inches='tight')
    print('plot saved in {0}'.format(energy_plot_filename))