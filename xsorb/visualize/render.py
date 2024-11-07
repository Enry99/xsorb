#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Enrico Pedretti

#This is a modified version of the render_image function from the xplot package,
#and will be sometimes updated/syncronized with it (they are basically the same).
#Now by putting this render_image as an isolated function it is very easy to
#keep it updated with the xplot package, since it is independent from the rest of the
#code and we just need to pass the right parameters to it.

#CHANGELOG:
# - 17 Oct 2024: taken from xplot, commit c0d4764.

"""
Module to render an image of an Atoms object using POVray or ASE renderer.

"""

from __future__ import annotations
import shutil
import os
from pathlib import Path

import numpy as np
from ase.io.pov import get_bondpairs
from ase.data.colors import jmol_colors as ATOM_COLORS
from ase import Atoms

from xsorb.ase_custom import AtomsCustom
from xsorb.ase_custom.io import write
import xsorb.ase_custom.povray  #monkey patching

#scaling factors for drawing atoms and bonds
ATOMIC_RADIUS_DEFAULT = 0.6
BOND_RADIUS_DEFAULT = 0.8
BOND_LINE_WIDTH_DEFAULT = 0.1


def _get_colorcoded_colors(atoms: Atoms, colorcode: str, ccrange : list | None = None) -> list:

    from matplotlib import colormaps, cm  #matplotlib==3.9.0
    from matplotlib.colors import Normalize

    if colorcode == 'forces':
        try:
            quantity = [np.linalg.norm(force) for force in atoms.get_forces()]
        except Exception as exc:
            raise ValueError("Forces are not present.") from exc
        cmap = colormaps.get_cmap('Blues')

    elif colorcode == 'magmoms':
        try:
            quantity = atoms.get_magnetic_moments()
        except Exception as exc:
            raise ValueError("Magnetic moments are not present.") from exc
        cmap = colormaps.get_cmap('coolwarm')

    elif colorcode == 'coordnum':
        from ase.neighborlist import NeighborList, natural_cutoffs
        cutoffs = natural_cutoffs(atoms, mult=1.2)
        nl = NeighborList(cutoffs, skin=0, self_interaction=False, bothways=True)
        nl.update(atoms)
        connmat = nl.get_connectivity_matrix()
        quantity = [connmat[idx].count_nonzero() for idx in range(len(atoms))]
        cmap = colormaps.get_cmap('viridis_r')
    else:
        raise ValueError("Invalid colorcode.")

    if ccrange is not None:
        vmin, vmax = ccrange[0], ccrange[1]
    else:
        vmin, vmax = min(quantity), max(quantity)
    norm = Normalize(vmin=vmin, vmax=vmax)
    scalar_map = cm.ScalarMappable(norm=norm, cmap=cmap)
    colors = [scalar_map.to_rgba(q)[:3] for q in quantity]

    return colors


def render_image(*,
                 atoms : Atoms | AtomsCustom,
                 outfile : str,
                 povray : bool = True,
                 width_res : int | None = 700,
                 rotations : str = '',
                 supercell : list | None = None,
                 wrap : bool = False,
                 transl_vector : list[float] | None = None,
                 depth_cueing : float | None = None,
                 range_cut : tuple | None = None,
                 cut_vacuum : bool = False,
                 colorcode : str | None = None,
                 ccrange : list | None = None,
                 arrows : str | None = None,
                 nobonds : bool = False,
                 mol_indices : list | None = None,
                 highlihgt_mol : bool = False,
                 custom_settings : dict | None = None):
    """
    Render an image of an Atoms object using POVray or ASE renderer.

    Args:
    - atoms: Atoms object to render.
    - outfile: path to the output file.
    - povray: if True, use POVray renderer (high quality, CPU intensive). If False, use ASE renderer
        (low quality, does not draw bonds).
    - width_res: width resolution of the output image.
    - rotations: string with the rotations to apply to the image.
    - supercell: list with the number of replicas in each direction.
        If mol_indices is provided, only the slab is replicated.
    - wrap: if True, wrap the atoms.
    - transl_vector: translation vector.
    - depth_cueing: intensity of depth cueing effect. If None, no depth cueing is applied.
    - range_cut: tuple with the range of z values to keep. If None, no range cut is applied.
    - cut_vacuum: if True, cut the vacuum in the z direction.
    - colorcode: if not None, color the atoms according to the specified quantity.
        Options are 'forces', 'magmoms' and 'coordnum'.
    - ccrange: list with the range of values to use for colorcoding.
        If None, the range is automatically set to the min and max of the quantity.
    - arrows: if not None, draw arrows for the specified quantity.
        Options are 'forces', 'magmoms' and 'coordnum'.
    - nobonds: if True, do not draw bonds.
    - mol_indices: list with the indices of the atoms to consider as the molecule
    - custom_settings: dictionary with custom settings, containing:
        'atomic_colors', 'atomic_radius', 'bond_radius', 'bond_line_width' and 'cell_line_width',
        and 'nontransparent_atoms'
    """

    atoms = atoms.copy() #do not modify the original object

    label = Path(outfile).stem

    if transl_vector is not None:
        atoms.translate(transl_vector)
        wrap = True

    if wrap:
        atoms.wrap()

    if supercell is not None:
        if mol_indices is not None:
            slab = Atoms([atom for i, atom in enumerate(atoms) if i not in mol_indices])
            slab *= supercell
            atoms = slab + Atoms([atom for i, atom in enumerate(atoms) if i in mol_indices])
            mol_indices = [i + len(slab) for i in mol_indices]
        else:
            atoms *= supercell

    if cut_vacuum:
        atoms.translate([0,0,atoms.positions[:,2].min()]) #shift to z=0
        atoms.cell[2,2] = atoms.positions[:,2].max() + 1
        atoms.pbc=[True,True,False] #to avoid periodic bonding in z direction

    if range_cut is not None:
        del atoms[[atom.index for atom in atoms if atom.z < range_cut[0] or atom.z > range_cut[1]]]
        atoms.translate([0,0,-range_cut[0]]) #shift the atoms to the origin of the new cell
        atoms.cell[2,2] = range_cut[1] - range_cut[0] #set the new cell height
        atoms.pbc=[True,True,False] #to avoid periodic bonding in z direction

    #set custom colors if present ################################################

    if custom_settings is not None:
        USER_COLORS       = custom_settings.get("atomic_colors", {})
        MOL_COLORS        = custom_settings.get("molecule_colors", {})
        ATOMIC_RADIUS     = custom_settings.get("atomic_radius", ATOMIC_RADIUS_DEFAULT)
        BOND_RADIUS       = custom_settings.get("bond_radius", BOND_RADIUS_DEFAULT)
        BOND_LINE_WIDTH   = custom_settings.get("bond_line_width", BOND_LINE_WIDTH_DEFAULT)
        CELLLINEWIDTH     = custom_settings.get("cell_line_width", 0)
    else:
        USER_COLORS  = {}
        MOL_COLORS   = {}
        ATOMIC_RADIUS     = ATOMIC_RADIUS_DEFAULT
        BOND_RADIUS       = BOND_RADIUS_DEFAULT
        BOND_LINE_WIDTH   = BOND_LINE_WIDTH_DEFAULT
        CELLLINEWIDTH     = 0


    if colorcode is None:
        #first, apply those of jmol
        colors = [ ATOM_COLORS[atom.number] for atom in atoms]

        #then, substitute user-defined colors
        try:
            species = atoms.custom_labels
        except AttributeError:
            species = atoms.get_chemical_symbols()

        for i, sp in enumerate(species):
            if mol_indices is not None and i in mol_indices:
                if sp in MOL_COLORS:
                    colors[i] = MOL_COLORS[sp]
                    continue
            if sp in USER_COLORS:
                colors[i] = USER_COLORS[sp]
    else:
        colors = _get_colorcoded_colors(atoms, colorcode, ccrange)


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

    if width_res is None:
        width_res = 700
        # I used 3000 in Xsorb paper. From 1500 is still quite good.
        # 2000 maybe best compromise (still very high res)

    if povray: #use POVray renderer (high quality, CPU intensive)
        config_copy = atoms.copy()
        #config_copy.set_pbc([0,0,0]) #to avoid drawing bonds with invisible replicas

        if mol_indices is not None and highlihgt_mol:
            textures = ['ase3' if i in mol_indices else 'pale' for i in range(len(atoms))]
        else:
            textures = None

        if custom_settings and custom_settings.get('nontransparent_atoms') is not None:
            transmittances = []
            textures = []
            trans_map = {True: 0.0, False: 0.8}
            texture_map = {True: 'ase3', False: 'pale'}
            for i in range(len(atoms)):
                nontrasp = i in custom_settings.get('nontransparent_atoms')
                transmittances.append(trans_map[nontrasp])
                textures.append(texture_map[nontrasp])
        else:
            transmittances = None

        povray_settings=dict(canvas_width=width_res,
                                celllinewidth=CELLLINEWIDTH,
                                transparent=False,
                                camera_type='orthographic',
                                camera_dist=2,
                                textures=textures,
                                transmittances=transmittances,
                                bondlinewidth=BOND_LINE_WIDTH,
                            )
        if not nobonds:
            povray_settings['bondatoms'] = get_bondpairs(config_copy, radius=BOND_RADIUS)
        if depth_cueing is not None:

            # Calculate height of ground fog by finding where the slab begins,
            # if we have a slab-molecule system

            if mol_indices is not None:
                mol_zs = atoms.positions[mol_indices][:,2]
                constant_fog_height = - (mol_zs.max() - mol_zs.min())
            else:
                # Try to guess by creating a histogram of z values, and considering that
                # the slab begins where the histogram is > 70% of the maximum density
                nbins = 5*int(atoms.positions[:,2].max() - atoms.positions[:,2].min()) +1 #5bin/Angstrom
                hist, bin_edges = np.histogram(atoms.positions[:,2], bins=nbins)
                threshold = 0.4 * hist.max()
                zmax_mol = atoms.positions[:,2].max()
                #Get highest bin where the histogram is above the threshold
                for i in range(len(hist)-1,-1,-1):
                    if hist[i] > threshold:
                        zmax_slab = bin_edges[i] + 1
                        break
                zmax_slab = min(zmax_mol, zmax_slab)
                constant_fog_height = - (zmax_mol - zmax_slab)

            povray_settings['depth_cueing'] = True
            povray_settings['cue_density'] = depth_cueing
            povray_settings['constant_fog_height'] = constant_fog_height


        #Do the actual rendering
        write(f'{label}.pov',
            atoms,
            format='pov',
            radii = ATOMIC_RADIUS,
            rotation=rotations,
            colors=colors,
            povray_settings=povray_settings,
            arrows_type = arrows
        ).render()
        os.remove(f'{label}.pov')
        os.remove(f'{label}.ini')

        if outfile != f'{label}.png': #move the output file to the desired location
            shutil.move(f'{label}.png', outfile)

    else: # use ASE renderer (low quality, does not draw bonds)
        write(outfile,
              atoms,
              format='png',
              radii = ATOMIC_RADIUS,
              rotation=rotations,
              colors=colors,
              maxwidth=width_res,
              scale=100)
