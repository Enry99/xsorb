#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue 28 Nov 2023

@author: Enrico Pedretti

Functions to render an image and to setup the rendering
"""


from ase.io.pov import get_bondpairs
from ase.io import read, write
from ase import Atoms
import json
import os, shutil
import numpy as np
from tqdm import tqdm
from render_atoms import ase_custom


def render_image(atoms : Atoms,
                 label : str,
                 outfile : str = None,
                 povray : bool = True,
                 width_res : int = 700,
                 rotations : str = '',
                 supercell : list = [1,1,1],
                 wrap : bool = False,
                 depth_cueing : float = None,
                 range_cut : tuple = None,
                 cut_vacuum : bool = False,
                 colorcode : str = None,
                 arrows : str = None,
                 nobonds : bool = False,
                 custom_settings = None):


    #scaling factors for drawing atoms and bonds
    ATOMIC_RADIUS_DEFAULT = 0.6
    BOND_RADIUS_DEFAULT = 0.8
    BOND_LINE_WIDTH_DEFAULT = 0.1

    if wrap:
        atoms.wrap()

    if supercell:
        atoms *= supercell

    if cut_vacuum:
        atoms.translate([0,0,-min([atom.z for atom in atoms])]) #shift to z=0
        atoms.cell[2,2] = np.max([atom.z for atom in atoms])
        atoms.pbc=[True,True,False] #to avoid periodic bonding in z direction

    if range_cut is not None:
        del atoms[[atom.index for atom in atoms if atom.z < range_cut[0] or atom.z > range_cut[1]]]
        atoms.translate([0,0,-range_cut[0]]) #shift the atoms to the origin of the new cell
        atoms.cell[2,2] = range_cut[1] - range_cut[0] #set the new cell height
        atoms.pbc=[True,True,False] #to avoid periodic bonding in z direction

    #set custom colors if present ################################################
    from ase.data.colors import jmol_colors
    ATOM_COLORS = jmol_colors.copy()

    if custom_settings is not None:
        USER_COLORS       = custom_settings["atomic_colors"] if "atomic_colors" in custom_settings else []
        ATOMIC_RADIUS     = custom_settings["atomic_radius"] if "atomic_radius" in custom_settings else ATOMIC_RADIUS_DEFAULT
        BOND_RADIUS       = custom_settings["bond_radius"] if "bond_radius" in custom_settings else BOND_RADIUS_DEFAULT
        BOND_LINE_WIDTH   = custom_settings["bond_line_width"] if "bond_line_width" in custom_settings else BOND_LINE_WIDTH_DEFAULT
        CELLLINEWIDTH     = custom_settings["cell_line_width"] if "cell_line_width" in custom_settings else 0
    else:
        USER_COLORS  = []
        ATOMIC_RADIUS     = ATOMIC_RADIUS_DEFAULT
        BOND_RADIUS       = BOND_RADIUS_DEFAULT
        BOND_LINE_WIDTH   = BOND_LINE_WIDTH_DEFAULT
        CELLLINEWIDTH     = 0


    if colorcode is None:
        #first, apply those of jmol
        colors = [ ATOM_COLORS[atom.number] for atom in atoms]

        #then, substitute user-defined colors
        try:
            species = atoms.get_custom_labels()
        except:
            species = atoms.get_chemical_symbols()

        for i, sp in enumerate(species):
            if sp in USER_COLORS:
                #print(sp)
                colors[i] = USER_COLORS[sp]
    else:
        from matplotlib import cm
        from matplotlib.colors import Normalize

        if colorcode == 'forces':
            try:
                quantity = [np.linalg.norm(force) for force in atoms.get_forces()]
            except:
                raise ValueError("Forces are not present.")
            cmap = cm.get_cmap('Blues')

        elif colorcode == 'magmoms':
            try:
                quantity = atoms.get_magnetic_moments()
            except:
                raise ValueError("Magnetic moments are not present.")
            cmap = cm.get_cmap('coolwarm')

        elif colorcode == 'coordnum':
            from ase.neighborlist import NeighborList, natural_cutoffs
            cutoffs = natural_cutoffs(atoms, mult=1.2)
            nl = NeighborList(cutoffs, skin=0, self_interaction=False, bothways=True)
            nl.update(atoms)
            connmat = nl.get_connectivity_matrix()
            quantity = [connmat[idx].count_nonzero() for idx in range(len(atoms))]
            cmap = cm.get_cmap('viridis_r')
        else:
            raise ValueError("Invalid colorcode.")

        norm = Normalize(vmin=min(quantity), vmax=max(quantity))
        scalar_map = cm.ScalarMappable(norm=norm, cmap=cmap)
        colors = [scalar_map.to_rgba(q)[:3] for q in quantity]


    ############################################################################
    # OLD depth cueing
    #fading color for lower layers in top view
    #if (depth_cueing is not None):
    #    zmax = max([atom.z for atom in atoms])
    #    zmin = min([atom.z for atom in atoms])
    #    delta = zmax - zmin
    #    if depth_cueing < 0:
    #        raise ValueError("depth_cueing_intensity must be >=0.")
    #    for atom in atoms:
    #        r,g,b = colors[atom.index] + (np.array([1,1,1]) - colors[atom.index])*(zmax - atom.z)/delta * depth_cueing
    #        if r>1: r=1
    #        if g>1: g=1
    #        if b>1: b=1
    #        colors[atom.index] = [r,g,b]
    ############################################################################


    if(povray): #use POVray renderer (high quality, CPU intensive)
        config_copy = atoms.copy()
        #config_copy.set_pbc([0,0,0]) #to avoid drawing bonds with invisible replicas

        povray_settings=dict(canvas_width=width_res,
                                celllinewidth=CELLLINEWIDTH,
                                transparent=False,
                                camera_type='orthographic',
                                camera_dist=2,
                                #textures = (['pale'] if depth_cueing else ['ase3']) * len(atoms),
                                bondlinewidth=BOND_LINE_WIDTH,
                            )
        if not nobonds:
            povray_settings['bondatoms'] = get_bondpairs(config_copy, radius=BOND_RADIUS)
        if depth_cueing is not None:
            povray_settings['depth_cueing'] = True
            povray_settings['cue_density'] = depth_cueing

        #print(atoms.positions[:,2])

        write('{0}.pov'.format(label),
            atoms,
            format='pov',
            radii = ATOMIC_RADIUS,
            rotation=rotations,
            colors=colors,
            povray_settings=povray_settings,
            arrows_type = arrows
        ).render()
        os.remove('{0}.pov'.format(label))
        os.remove('{0}.ini'.format(label))

        if outfile:
            #we have to do this here since we cannot just set the output path in write_pov,
            #since the ini and pov files must be in the same folder, without absolute paths
            shutil.move('{0}.png'.format(label), outfile)

    else: # use ASE renderer (low quality, does not draw bonds)
        write(outfile if outfile else label + '.png',
              atoms,
              format='png',
              radii = ATOMIC_RADIUS,
              rotation=rotations,
              colors=colors,
              maxwidth=width_res,
              scale=100)




def start_rendering(filename : str,
                    index : str = '-1',
                    movie : bool = False,
                    framerate : int = 10,
                    custom_settings : dict = None,
                    **kwargs
                    ):


    if index == '-1' and movie: #if we want to render a movie, we need to read the whole trajectory
        index = ':'
    atoms = read(filename, index=index) #read file (any format supported by ASE)
    label = os.path.splitext(os.path.basename(filename))[0]  #get filename without extension, and without path before possible /
    print('File was read successfully.')

    if type(atoms) is list:

        print(f'Rendering {len(atoms)} images...')

        os.makedirs('rendered_frames', exist_ok=True)
        main_dir = os.getcwd()
        os.chdir('rendered_frames')

        for i, atom in enumerate(tqdm(atoms)):
            render_image(atom,
                         label + '_{:05d}'.format(i),
                         custom_settings=custom_settings,
                         **kwargs)
        print('Rendering complete.')

        os.chdir(main_dir)

        if movie:
            print('Generating movie...')
            os.system(f'ffmpeg -framerate {framerate} -i rendered_frames/{label}_%05d.png  -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2"  -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p {label}.mp4')
            print('Movie generated.')
    else:
        print('Rendering image...')
        render_image(atoms,
                     label=label,
                     custom_settings=custom_settings,
                     **kwargs)
        print('Rendering complete.')

    print('Job done.')