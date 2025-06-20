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
import subprocess
import logging
from dataclasses import asdict

import xsorb.structures.slab
from xsorb.ase_custom.atoms import AtomsCustom
from xsorb.settings import Settings
from xsorb.io.database import Database
from xsorb.ase_custom.io import ase_custom_read as read
from xsorb.visualize.render import render_image
from xsorb.visualize.plot import plot_overview_grid
from xsorb.visualize.utils import get_centered_mol_and_slab
from xsorb.io.utils import progressbar
from xsorb.adsorptiondata.adsorptioncalculation import ALLOWED_STATUSES
from xsorb.adsorptiondata import AdsorptionStructure, AdsorptionCalculation
from xsorb.visualize.settings import CustomSettings


stars_map = {"completed": "",
             "incomplete": "*",
             "scf_nonconverged": "**"}
symbols_map = {"completed": {"symbol": ".", "color": "white"}, # not to be used
               "incomplete": {"symbol": "x", "color": "black"},
               "scf_nonconverged": {"symbol": "^", "color": "red"}}
assert stars_map.keys() == symbols_map.keys() and set(stars_map.keys()) == set(ALLOWED_STATUSES) \
    , "stars_map and symbols_map must have the same keys as ALLOWED_STATUSES"


def plot_adsorption_sites(all_sites : bool = False):
    '''
    Plot an image of the surface with the adsorption sites.

    Args:
    - all_sites: if True, plot all sites ignoring symm_reduce and selected_sites
    '''

    settings = Settings()

    slab = read(settings.input.slab_filename)

    slab = xsorb.structures.slab.Slab(slab=slab,
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
                framerate : int = 10,
                **kwargs):
    '''
    Plot images of the configurations

    Args:
    - calc_type: 'initial','screening','relax','mlopt'
    - calc_id: index of the calculation to plot. If None, plot all
    - movie: if True, generate a movie from the images
    - framerate: framerate for the movie, Default is 10.

    kwargs are those for xsorb.visualize.render.render_image
    '''

    kwargs.pop('traceback', None)
    kwargs.pop('command', None)
    kwargs.pop('func', None)
    framerate = kwargs.pop('framerate', None)
    rotations = kwargs.pop('rotation', None)
    center_mol = kwargs.pop('center_molecule', None)

    if calc_type not in ('initial','screening', 'relax', 'mlopt'):
        raise RuntimeError(f"Wrong '{calc_type}', expected 'screening', 'relax' or 'mlopt'")

    if calc_type == 'initial':
        rows = Database.get_structures(calc_ids=calc_id)
    else:
        rows = Database.get_calculations(calc_type, calc_ids=calc_id)

    if not rows:
        print("No images to be generated.")
        return

    custom_settings = CustomSettings()

    #get it here, so that we can decide to apply it or not
    #depending on the rotation.
    depth_cueing = kwargs.pop('depth_cueing', None)


    main_dir = os.getcwd()
    figures_dir = Path(f"{calc_type}_images/").absolute().as_posix()
    os.makedirs(figures_dir, exist_ok=True)
    os.chdir(figures_dir)

    outfiles = []
    energies = []
    stars = [] #for marking non-converged calculations
    for row in progressbar(rows, 'Rendering:'):

        atoms = AtomsCustom(row.toatoms())
        mol_indices = row.data.AdsorptionStructure.mol_indices

        if center_mol:
            #use the same translation for all frames in the trajectory, to avoid jumps
            _, _, transl_vector = get_centered_mol_and_slab(atoms, mol_indices)
        else:
            transl_vector = None

        if not rotations: #lateral, top,
            rot_list = ['', '-5z,-85x']
            rotations_labels = ['top', 'lateral']
            depth_cueings = [depth_cueing, None]
        else: #use specified rotation
            if 'front' in rotations:
                rotations = '-90x'
                rotation_label = 'front'
            else:
                rotation_label = rotations.replace(',','_')
            rot_list = [rotations]
            rotations_labels = [rotation_label]
            depth_cueings = [depth_cueing]

        for rot, rot_label, dc in zip(rot_list, rotations_labels, depth_cueings):

            file_label = f'{calc_type}_{row.calc_id}_{rot_label}'

            #render trajectory for each config
            if movie and (rot_label == rotations_labels[1] if not rotations else True):

                print('Generating frames for traj...')

                os.makedirs(f'rendered_frames_{row.calc_id}', exist_ok=True)
                os.chdir(f'rendered_frames_{row.calc_id}')

                adscalc = AdsorptionCalculation.fromdict(row.data["AdsorptionCalculation"])
                traj = adscalc.calc_results.trajectory
                for i, frame in enumerate(progressbar(traj, 'Rendering:')):
                    render_image(atoms=frame,
                                outfile=f'{file_label}_{i:05d}.png',
                                rotations=rot,
                                transl_vector=transl_vector,
                                depth_cueing=dc,
                                mol_indices=mol_indices,
                                custom_settings=custom_settings,
                                **kwargs)
                os.chdir(figures_dir)

                print('Frames generated. Generating movie...')

                success = False
                # first, try to use ffmpeg:
                ffmpeg_cmd = f'ffmpeg -y -framerate {framerate} '\
                    f'-i rendered_frames/{file_label}_%05d.png '\
                    '-vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" '\
                    f'-c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p '\
                    f'{file_label}.mp4'
                try:
                    ret = subprocess.run([ffmpeg_cmd], check=True, capture_output=True, shell=True)
                    success = ret.returncode == 0
                except (subprocess.CalledProcessError, FileNotFoundError) as ffmpeg_error:
                    logging.error('ffmpeg failed: %s. trying to use imagemagick...', ffmpeg_error)
                    # if ffmpeg fails, try imagemagick:
                    try:
                        convert_cmd = f'convert -delay {1000 // framerate} '\
                            f'-loop 0 rendered_frames/{file_label}_*.png {file_label}.gif'
                        ret = subprocess.run([convert_cmd], check=True, capture_output=True, shell=True)
                        success = ret.returncode == 0
                    except (subprocess.CalledProcessError, FileNotFoundError) as imagemagick_error:
                        logging.error('Imagemagick also failed: %s', imagemagick_error)

                if success:
                    logging.info('Movie generated.')
                else:
                    logging.error('Error generating movie, '
                            'however the frames are still present in the rendered_frames folder.')

            else: #single image
                render_image(atoms=atoms,
                         outfile=f'{file_label}.png',
                         rotations=rot,
                         transl_vector=transl_vector,
                         custom_settings=custom_settings,
                         depth_cueing=dc,
                         mol_indices=mol_indices,
                         **kwargs)
                if rot_label == rotations_labels[0]:
                    outfiles.append(f'{file_label}.png')
                    energies.append(row.get('adsorption_energy'))
                    stars.append('**' if row.get('scf_nonconverged') else \
                             '*' if row.get('status') != 'completed' else '')

    if calc_id is None and not movie:   #plot grid
        plot_overview_grid(calc_type=calc_type,
                           outfiles=outfiles,
                           rot_label=rotations_labels[0],
                           calc_indices=[row.calc_id for row in rows],
                           energies=energies,
                           stars=stars)
    #rotations_labels[0]: use top if not customrot, else use customrot for the grid

    os.chdir(main_dir)

    print(f'All images saved in {figures_dir}.')


def view_config(calc_type : str, calc_id : int):
    '''
    View the selected config with ASE GUI

    Args:
    - calc_type: 'initial','screening','relax','mlopt'
    - calc_id: index of the calculation to plot.
    '''

    from ase.visualize import view # pylint: disable=import-outside-toplevel

    if calc_type not in ('initial','screening', 'relax', 'mlopt'):
        raise RuntimeError(f"Wrong '{calc_type}', expected '"\
                           "'screening', 'relax', 'mlopt' or 'initial'")

    if calc_type == 'initial':
        row = Database.get_structures(calc_ids=calc_id)[0]
        ads_structure = AdsorptionStructure.fromdict(row.data["AdsorptionStructure"])
        atoms_or_traj = ads_structure.atoms
    else:
        row = Database.get_calculations(calc_type=calc_type, calc_ids=calc_id)[0]
        calc_results = AdsorptionCalculation.fromdict(row.data["AdsorptionCalculation"])
        atoms_or_traj = calc_results.calc_results.trajectory

    view(atoms_or_traj)


def plot_energy_evolution(calc_type : str):
    '''
    Plot the energy evolution during optimization for all the configurations
    that have at least one step.

    Args:
    - calc_type: 'screening','relax','mlopt'
    '''

    from matplotlib import pyplot as plt
    import numpy as np

    rows = Database.get_calculations(calc_type, selection='adsorption_energy')

    for row in rows:

        label = f'{row.calc_id}: {row.adsorption_energy:.2f}{stars_map[row.status]} eV'

        energy_array = row.data.AdsorptionCalculation.calc_results.adsorption_energy_evol
        plt.plot(energy_array, '-', label=label)

        if row.status != 'completed':
            plt.plot(len(energy_array)-1, energy_array[-1],
                     symbols_map[row.status]["sym"],
                     color=symbols_map[row.status]["color"])


    plt.title(f'Energy evolution during {calc_type}')
    plt.xlabel('step')
    plt.ylabel('energy (eV)')
    plt.grid(linestyle='dotted')
    plt.legend(title="Config, energy",
               ncols=np.ceil(len(rows)/10),
               prop={'size': 6  if calc_type in ('screening', 'mlopt') else 8})
    energy_plot_filename = f'{calc_type}_energies.png'
    plt.savefig(energy_plot_filename, dpi=300, bbox_inches='tight')
    print(f'Plot saved in {energy_plot_filename}')
