#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Enrico Pedretti

"""
Collection of functions for generating images that
are called from the command line.

"""

from __future__ import annotations
from pathlib import Path
import os
from dataclasses import asdict

from ase.visualize import view

import xsorb.structures.slab
from xsorb.io.settings import Settings
from xsorb.io.database import Database
from xsorb.visualize.render import render_image
from xsorb.visualize.plot import plot_overview_grid
from xsorb.visualize.utils import get_centered_mol_and_slab, read_custom_colors
from xsorb.io.utils import progressbar


def plot_adsorption_sites(all_sites : bool = False):
    '''
    Plot an image of the surface with the adsorption sites.

    Args:
    - all_sites: if True, plot all sites ignoring symm_reduce and selected_sites
    '''

    settings = Settings()

    slab = xsorb.structures.slab.Slab(slab_filename=settings.input.slab_filename,
            surface_thickness=settings.structure.adsorption_sites.surface_thickness,
            layers_threshold=settings.structure.constraints.layers_height,
            sort_atoms_by_z=settings.structure.misc.sort_atoms_by_z,
            translate_slab_from_below_cell_bottom=settings.structure.misc.translate_slab)

    #Find adsorption sites and labels (site type and x,y coords.)
    sites_settings = settings.structure.adsorption_sites

    if all_sites:
        if sites_settings.high_symmetry_params:
            sites_settings.high_symmetry_params.symm_reduce = 0.0
        sites_settings.selected_sites = None

    mode_params = {
        'high_symmetry': sites_settings.high_symmetry_params,
        'coord_number': sites_settings.coord_number_params,}
    mode = sites_settings.mode
    if mode not in mode_params:
        raise ValueError(f"mode must be one of {mode_params.keys()}")

    slab.find_adsorption_sites(
        mode=mode,
        **asdict(mode_params[mode]),
        selected_sites=sites_settings.selected_sites,
        save_image=True,
        verbose=True)


def plot_images(calc_type : str,
                calc_id : int | None = None,
                movie: bool = False,
                **kwargs):
    '''
    Plot images of the configurations

    Args:
    - calc_type: 'initial','screening','relax','ml_opt'
    - calc_id: index of the calculation to plot. If None, plot all
    - movie: if True, generate a movie from the images

    kwargs are those for xsorb.visualize.render_image
    '''

    if calc_type not in ('initial','screening', 'relax', 'ml_opt'):
        raise RuntimeError(f"Wrong '{calc_type}', expected 'screening', 'relax' or 'ml_opt'")

    if calc_type == 'initial':
        rows = Database.get_structures(calc_ids=calc_id)
    else:
        rows = Database.get_calculations(calc_type, calc_ids=calc_id)

    if not rows:
        print("No images to be generated.")
        return

    custom_colors = read_custom_colors()
    if custom_colors.get('povray_old_style'):
        os.environ['POVRAY_OLD_STYLE'] = '1'

    #get it here, so that we can decide to apply it or not
    #depending on the rotation.
    depth_cueing = kwargs.pop('depth_cueing')


    main_dir = os.getcwd()
    figures_dir = Path(f"{calc_type}_images/").absolute().as_posix()
    os.makedirs(figures_dir, exist_ok=True)
    os.chdir(figures_dir)

    outfiles = []
    energies = []
    stars = [] #for marking non-converged calculations
    for row in progressbar(rows, 'Rendering:'):

        atoms = row.to_atoms()
        mol_indices = row.data.get('mol_indices')

        if kwargs.get('center_mol'):
            #use the same translation for all frames in the trajectory, to avoid jumps
            _, _, transl_vector = get_centered_mol_and_slab(atoms, mol_indices)
        else:
            transl_vector = None

        if not kwargs.get('rotations'): #lateral, top,
            rot_list = ['-5z,-85x', '']
            rotations_labels = ['lateral', 'top', ]
            depth_cueings = [None, depth_cueing]
        else: #use specified rotation
            rot_list = [kwargs.get('rotations')]
            rotations_labels = [kwargs.get('rotations').replace(',','_')]
            depth_cueings = [depth_cueing]

        for rot, rot_label, dc in zip(rot_list, rotations_labels, depth_cueings):

            file_label = f'{calc_type}_{row.calc_id}_{rot_label}'

            #render trajectory for each config
            if movie and rot_label == rotations_labels[0]:

                print('Generating frames for traj...')

                os.makedirs(f'rendered_frames_{row.calc_id}', exist_ok=True)
                os.chdir(f'rendered_frames_{row.calc_id}')
                for i, frame in enumerate(progressbar(row.data.trajectory,
                                                        'Rendering:')):
                    render_image(atoms=frame,
                                outfile=f'{file_label}_{i:05d}.png',
                                rotations=rot,
                                transl_vector=transl_vector,
                                custom_settings=custom_colors,
                                depth_cueing=dc,
                                **kwargs)
                os.chdir(figures_dir)

                print('Frames generated. Generating movie...')

                success = os.system(f'ffmpeg -framerate {kwargs.get("framerate")} '\
                            f'-i rendered_frames/{file_label}_%05d.png  '\
                            '-vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" '\
                            f'-c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p '\
                            f'{file_label}.mp4')

                if success:
                    print('Movie generated.')
                else:
                    print('Error generating movie, '\
                          'however the frames are still present in the  rendered_frames folder.')

            else: #single image
                render_image(atoms=atoms,
                         outfile=f'{file_label}.png',
                         rotations=rot,
                         transl_vector=transl_vector,
                         custom_settings=custom_colors,
                         depth_cueing=dc,
                         **kwargs)
                outfiles.append(f'{file_label}.png')
                energies.append(row.get('adsorption_energy'))
                stars.append('**' if row.get('scf_nonconverged') else \
                             '*' if row.get('status') != 'completed' else '')

    if calc_id is None and not movie:   #plot grid
        plot_overview_grid(calc_type=calc_type,
                           outfiles=outfiles,
                           calc_indices=[row.calc_id for row in rows],
                           energies=energies,
                           stars=stars)
    #rotations_labels[0]: use top if not customrot, else use customrot for the grid

    os.chdir(main_dir)

    print(f'All images saved in {figures_dir}.')


#TODO: make it read trajectory
def view_config(calc_type : str, calc_id : int):
    '''
    View the selected config with ASE GUI

    Args:
    - calc_type: 'initial','screening','relax','ml_opt'
    - calc_id: index of the calculation to plot.
    '''

    if calc_type not in ('initial','screening', 'relax', 'ml_opt'):
        raise RuntimeError(f"Wrong '{calc_type}', expected 'screening', 'relax' or 'ml_opt'")

    if calc_type == 'initial':
        atoms = Database.get_structure(calc_id)[0].to_atoms()
    else:
        atoms = Database.get_calculation(calc_type, calc_id)[0].to_atoms()

    view(atoms)


def plot_energy_evolution(calc_type : str):
    '''
    Plot the energy evolution during optimization for all the configurations

    Args:
    - calc_type: 'screening','relax','ml_opt'
    '''

    from matplotlib import pyplot as plt
    import numpy as np

    rows = Database.get_calculations(calc_type, selection='adsorption_energy')

    for row in rows:

        stars = '**' if row.get('scf_nonconverged') else \
            '*' if row.get('status') != 'completed' else ''

        label = f'{row.calc_id}: {row.get("adsorption_energy"):.2f}{stars} eV'

        energy_array = row.data.adsorption_energy_evolution
        plt.plot(energy_array, '-', label=label)

        if '*' in stars:
            symbol = 'x' if stars == '*' else '^'
            color = 'black' if stars == '*' else 'red'
            plt.plot(len(energy_array)-1, energy_array[-1], symbol, color=color)


    plt.title('Energy evolution during optimization')
    plt.xlabel('step')
    plt.ylabel('energy (eV)')
    plt.grid(linestyle='dotted')
    plt.legend(title="Config, energy",
               ncols=np.ceil(len(rows)/10),
               prop={'size': 6  if calc_type in ('screening', 'ml_opt') else 8})
    energy_plot_filename = f'{calc_type}_energies.png'
    plt.savefig(energy_plot_filename, dpi=300, bbox_inches='tight')
    print(f'Plot saved in {energy_plot_filename}')
