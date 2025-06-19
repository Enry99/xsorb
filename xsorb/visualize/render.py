#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Enrico Pedretti

#This is a modified version of the render_image function from the atomsplot package,
#and will be sometimes updated/syncronized with it (they are basically the same).
#Now by putting this render_image as an isolated function it is very easy to
#keep it updated with the atomsplot package, since it is independent from the rest of the
#code and we just need to pass the right parameters to it.

#CHANGELOG:
# - 19 Jun 2025: sync with atomsplot v1.0

"""
Module to render an image of an Atoms object using POVray or ASE renderer.

"""

from __future__ import annotations
import logging
from typing import Optional, TYPE_CHECKING
import shutil
import os
from pathlib import Path

import numpy as np
from ase.io.pov import get_bondpairs, set_high_bondorder_pairs
from ase.data import covalent_radii

from ase.io import write
from ase.io.utils import PlottingVariables
from ase.io.pov import POVRAY, POVRAYIsosurface
from ase.geometry.geometry import get_layers
from ase.build.tools import sort

from xsorb.visualize.settings import CustomSettings
from xsorb.ase_custom import AtomsCustom # monkey patch for ase.utils.PlottingVariables arrows_type. pylint: disable=unused-import
import xsorb.ase_custom.povray # monkey patch for povray. pylint: disable=unused-import

if TYPE_CHECKING:
    from ase import Atoms
    from xsorb.ase_custom import AtomsCustom



def _get_colorcoded_colors(atoms: Atoms, quantity: str, ccrange : list | None = None) -> list:
    """
    Get colors for atoms based on the specified quantity.

    Parameters
    ----------
    atoms : Atoms
        Atoms object containing the atomic structure.
    quantity : str
        The property to color the atoms by. Options are 'forces', 'magmoms', 'coordnum'.
    ccrange : list, optional
        Range of values for color coding. If None, the range is automatically set to the
        min and max of the quantity.

    Returns
    -------
    colors : list
        List of RGB colors for each atom, normalized to the specified range.
    """

    # import here to reduce loading time when not needed
    from matplotlib import colormaps, cm  #matplotlib==3.9.0 # pylint: disable=import-outside-toplevel
    from matplotlib.colors import Normalize # pylint: disable=import-outside-toplevel

    if quantity == 'forces':
        try:
            values = np.linalg.norm(atoms.get_forces(), axis=1)
        except Exception as exc:
            raise ValueError("Forces are not present.") from exc
        cmap = colormaps.get_cmap('Blues')

    elif quantity == 'magmoms':
        try:
            values = atoms.get_magnetic_moments()
        except Exception as exc:
            raise ValueError("Magnetic moments are not present.") from exc
        cmap = colormaps.get_cmap('coolwarm')

    elif quantity == 'coordnum':
        from ase.neighborlist import NeighborList, natural_cutoffs # pylint: disable=import-outside-toplevel
        cutoffs = natural_cutoffs(atoms, mult=1.3)
        nl = NeighborList(cutoffs, skin=0, self_interaction=False, bothways=True)
        nl.update(atoms)
        connmat = nl.get_connectivity_matrix()
        values = [connmat[idx].count_nonzero() for idx in range(len(atoms))]
        cmap = colormaps.get_cmap('viridis_r')
    else:
        raise ValueError("Invalid quantity for colorcoding.")

    if ccrange is not None:
        vmin, vmax = ccrange[0], ccrange[1]
    else:
        vmin, vmax = min(values), max(values)
    norm = Normalize(vmin=vmin, vmax=vmax)
    scalar_map = cm.ScalarMappable(norm=norm, cmap=cmap)
    colors = [scalar_map.to_rgba(q)[:3] for q in values]

    return colors


def _get_arrows(atoms: Atoms, quantity: str, rotation, scaling : float) -> np.ndarray:

    if quantity == 'forces':
        arrows = atoms.get_forces()
    elif quantity == 'magmoms':
        arrows = np.array([[0,0, magmom] for magmom in atoms.get_magnetic_moments()])
    else:
        raise ValueError(f'Unknown arrows type: {quantity}')

    #normalize arrows
    maxlength = np.linalg.norm(arrows, axis=1).max()
    arrows = arrows / maxlength

    arrows = np.dot(arrows, rotation) * scaling

    return arrows


def _calculate_ground_fog_height(atoms: Atoms, mol_indices: list[float] | None = None) -> float:
    """
    Calculate the height of ground fog for visualization purposes.
    This function determines the fog height by analyzing the structure of a system,
    particularly distinguishing between slab and molecular components.
    The fog height is used for depth cueing in visualizations, providing a sense of depth.

    Parameters
    ----------
    atoms : Atoms
    mol_indices : list[float] | None, optional
        List of indices corresponding to molecular atoms. If provided, the fog
        height is calculated based on the z-coordinate range of these atoms.
        If None, the function attempts to distinguish between slab and molecular
        regions automatically.

    Returns
    -------
    float
        The calculated ground fog height as a negative value representing the
        vertical extent from the top of the slab to the top of the molecules.

    Notes
    -----
    When mol_indices is None:
        - The function sorts atoms by z-coordinate and uses layer analysis to
          identify slab regions.
        - Layers are identified using Miller indices (0,0,1) with a tolerance of 0.3.
        - Volume-based thresholding (80% of max occupied volume) is used to
          distinguish slab layers from molecular regions.
        - Occupied volume is estimated using covalent radii assuming spherical atoms.
    """

    if mol_indices is not None:
        mol_zs = atoms.positions[mol_indices][:,2]
        constant_fog_height = - (mol_zs.max() - mol_zs.min())
    else:
        #this should give much more intense peaks for the slab
        sorted_atoms = sort(atoms, tags=atoms.positions[:,2])
        zmax_mol = sorted_atoms.positions[-1,2]
        layer_indicization, distances = get_layers(atoms=sorted_atoms, miller=(0,0,1),
                                        tolerance=0.3)

        #use occupied volume rather than number of atoms, to avoid
        #planar molecules such as benzene to be considered as a slab layer
        bins_heights = [0] * len(distances)
        for i, idx in enumerate(layer_indicization):
            bins_heights[idx] += 4/3*np.pi*(covalent_radii[sorted_atoms[i].number])**3

        # XXX: this method does not work for the example depth_cueing_mols_auto
        threshold = 0.8 * max(bins_heights)
        zmax_slab = 0
        for i in range(len(distances)-1,-1,-1):
            if bins_heights[i] > threshold:
                zmax_slab = distances[i] +1
                break

        if zmax_slab < zmax_mol - 5:
            logging.warning(
                "Not possible to determine slab height, resorting to default fog "
            )
            constant_fog_height = -2
        else:
            zmax_slab = min(zmax_mol, zmax_slab) #avoid negative heights
            constant_fog_height = - (zmax_mol - zmax_slab)

    return constant_fog_height


def _calculate_bondorder_pairs(atoms: Atoms, mol_indices : list[int] | None = None) -> dict:
    '''
    Calculate pairs of atoms with high bond order using Pymatgen's CovalentBondNN.
    '''

    from pymatgen.analysis.local_env import CovalentBondNN # pylint: disable=import-outside-toplevel
    from pymatgen.io.ase import AseAtomsAdaptor # pylint: disable=import-outside-toplevel

    covalent_bond_nn = CovalentBondNN()
    pymat_mol = AseAtomsAdaptor.get_molecule(atoms)

    if mol_indices is None:
        mol_indices = list(range(len(pymat_mol)))

    high_bondorder_pairs = {}
    for atom_id in mol_indices:
        neighbors = covalent_bond_nn.get_nn_info(pymat_mol, atom_id)

        for neighbor in neighbors:
            bondorder = round(neighbor['weight'])
            # avoid duplicate pairs
            if  bondorder >= 2 and (neighbor['site_index'], atom_id) \
                not in high_bondorder_pairs:
                high_bondorder_pairs[(atom_id, neighbor['site_index'])] = \
                    ((0, 0, 0), bondorder, (0.2, 0.2, 0))

    return high_bondorder_pairs


def render_image(atoms: 'Atoms | AtomsCustom',
                outfile: str,
                custom_settings: CustomSettings,
                rotations: str = '',
                supercell: Optional[list] = None,
                wrap: bool = False,
                range_cut: Optional[tuple] = None,
                cut_vacuum: bool = False,
                bonds: str = 'single',
                hide_cell: bool = False,
                depth_cueing: Optional[float] = None,
                highlight_mol: bool = False,
                colorcode: Optional[str] = None,
                ccrange: Optional[list] = None,
                arrows: Optional[str] = None,
                arrows_scale: float = 1.0,
                chg_grid: Optional[np.ndarray] = None,
                chg_iso_threshold: Optional[float] = None,
                width_res: Optional[int] = 700,
                povray: bool = True,
                transl_vector: Optional[list[float]] = None,
                mol_indices: Optional[list] = None,
                fixed_bounds : bool = False):

    """
    Render an image of an Atoms object using POVray or ASE renderer.

    Parameters
    ----------
    atoms : Atoms | AtomsCustom
        Atoms object to render.
    outfile : str
        Path to the output file.
    custom_settings : CustomSettings
        Custom settings for rendering, including colors, radii, and other parameters.
    rotations : str, optional
        String with the rotations to apply to the image. Default is ''.
    supercell : list | None, optional
        List with the number of replicas in each direction.
        If mol_indices is provided, only the slab is replicated. Default is None.
    wrap : bool, optional
        If True, wrap the atoms. Default is False.
    range_cut : tuple | None, optional
        Tuple with the range of z values to keep. If None, no range cut is applied. Default is None.
    cut_vacuum : bool, optional
        If True, cut the vacuum in the z direction. Default is False.
    bonds : str, optional
        Type of bonds to draw. Options are 'none', 'single' (default), 'multiple'.
    hide_cell : bool, optional
        If True, hide the cell box. Default is False.
    depth_cueing : float | None, optional
        Intensity of depth cueing effect. If None, no depth cueing is applied. Default is None.
    highlihgt_mol : bool, optional
        If True, highlight molecular atoms with a different color. Default is False.
    colorcode : str | None, optional
        If not None, color the atoms according to the specified quantity.
        Options are 'forces', 'magmoms' and 'coordnum'. Default is None.
    ccrange : list | None, optional
        List with the range of values to use for colorcoding.
        If None, the range is automatically set to the min and max of the quantity. Default is None.
    arrows : str | None, optional
        If not None, draw arrows for the specified quantity.
        Options are 'forces', 'magmoms' and 'coordnum'. Default is None.
    arrows_scale : float, optional
        Scale factor for the arrows. Default is 1.0 (no scaling).
    chg_grid : np.ndarray | None, optional
        Charge density grid to use for isosurface rendering.
        If None, no isosurface is rendered. Default is None.
    chg_iso_threshold : float | None, optional
        Iso-surface threshold for the charge density.
        If None, VESTA default is used (mean(|rho|) + 2 * std(|rho|)). Default is None.
    width_res : int | None, optional
        Width resolution of the output image. Default is 700.
    povray : bool, optional
        If True, use POVray renderer (high quality, CPU intensive). If False, use ASE renderer
        (low quality, does not draw bonds). Default is True.
    transl_vector : list[float] | None, optional
        Translation vector for the molecule. Default is None.
    mol_indices : list | None, optional
        List with the indices of the atoms to consider as the molecule. Default is None.
    """

    calc = atoms.calc
    atoms = atoms.copy() #do not modify the original object
    atoms.calc = calc #keep the calculator, if present

    label = Path(outfile).stem

    if transl_vector is not None:
        atoms.translate(transl_vector)
        wrap = True

    if wrap:
        atoms.wrap(pretty_translation=True) #wrap the atoms to the unit cell

    if mol_indices is None and custom_settings.mol_indices is not None:
        mol_indices = custom_settings.mol_indices

    if supercell is not None:
        if mol_indices is not None:
            n_repeats = np.array(supercell) - 1
            # get new mol_indces after supercell expansion
            mol_indices = [i + j * len(atoms) for i in mol_indices for j in range(np.prod(supercell))]

        #     slab_indices = [i for i in range(len(atoms)) if i not in mol_indices]
        #     slab = atoms[slab_indices]
        #     slab *= supercell
        #     atoms = slab + atoms[mol_indices]
        #     mol_indices = [i + len(slab) for i in range(len(mol_indices))]
        # else:
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

    #set custom colors if present ###############################################


    if colorcode is None:
        #first, apply those of jmol
        colors = [ custom_settings.color_scheme[atom.number] for atom in atoms]

        #then, substitute user-defined colors
        if isinstance(atoms, AtomsCustom):
            species = atoms.custom_labels
        else:
            species = atoms.get_chemical_symbols()

        for i, sp in enumerate(species):
            if mol_indices is not None and i in mol_indices:
                if sp in custom_settings.molecule_colors:
                    colors[i] = custom_settings.molecule_colors[sp]
                    continue
            if sp in custom_settings.atomic_colors:
                colors[i] = custom_settings.atomic_colors[sp]
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
        # I used 3000 in Xsorb paper. > 1500 is still very good.

    if povray: #use POVray renderer (high quality, CPU intensive)
        pvars = PlottingVariables(atoms,
            scale=1,
            radii=custom_settings.atomic_radius,
            rotation=rotations,
            colors=colors,
            auto_bbox_size=1.2 if fixed_bounds else 1.05, #auto_bbox_size is used to set the size of the bounding box
            show_unit_cell=3 if fixed_bounds else 2,  #IMPORTANT: keep the blank space around the cell fixed in trajs
        )

        if mol_indices is not None and highlight_mol:
            textures = ['ase3' if i in mol_indices else 'pale' for i in range(len(atoms))]
        else:
            textures = None

        if custom_settings.nontransparent_atoms:
            transmittances = []
            textures = []
            trans_map = {True: 0.0, False: 0.8}
            texture_map = {True: 'ase3', False: 'pale'}
            for i in range(len(atoms)):
                nontrasp = i in custom_settings.nontransparent_atoms
                transmittances.append(trans_map[nontrasp])
                textures.append(texture_map[nontrasp])
        else:
            transmittances = None

        if 'x' not in rotations and 'y' not in rotations:
            dz = atoms.cell[2,2] - atoms.positions[:,2].max() + 0.1
        else:
            dz = 0
        camera_dist = max(2, dz)


        povray_settings=dict(
            canvas_width=width_res,
            celllinewidth=custom_settings.cell_line_width if not hide_cell else 0,
            transparent=False,
            camera_type='orthographic',
            camera_dist=camera_dist,
            textures=textures,
            transmittances=transmittances,
            bondlinewidth=custom_settings.bond_line_width,
            arrows = _get_arrows(atoms, arrows, pvars.rotation, arrows_scale)
            if arrows is not None else None
        )

        if bonds == 'none' or custom_settings.nontransparent_atoms:
            # with transparency, bonds are very ugly
            bondatoms = None
        else:
            bondatoms = get_bondpairs(atoms, radius=custom_settings.bond_radius)
            if bonds == 'multiple':
                high_bondorder_pairs = _calculate_bondorder_pairs(atoms)
                bondatoms = set_high_bondorder_pairs(bondatoms, high_bondorder_pairs)
            povray_settings['bondatoms'] = bondatoms


        if depth_cueing is not None:
            constant_fog_height = _calculate_ground_fog_height(atoms) #, mol_indices)

            povray_settings['depth_cueing'] = True
            povray_settings['cue_density'] = depth_cueing
            povray_settings['constant_fog_height'] = constant_fog_height



        pov_obj = POVRAY.from_PlottingVariables(pvars, **povray_settings)


        if chg_grid is not None:
            if chg_iso_threshold is None:
                # VESTA default isosurface: mean(|rho|) + 2 * std(|rho|)
                chg_iso_threshold = np.mean(np.abs(chg_grid)) + 2 * np.std(np.abs(chg_grid))

            iso_positive = POVRAYIsosurface.from_POVRAY(
                povray=pov_obj,
                density_grid=chg_grid,
                cut_off=chg_iso_threshold,
                color=(0.80, 0.80, 0.0, 0.3))

            iso_negative = POVRAYIsosurface.from_POVRAY(
                povray=pov_obj,
                density_grid=chg_grid,
                cut_off=-chg_iso_threshold,
                color=(0.00, 0.80, 0.80, 0.3))

            pov_obj.isosurfaces = [iso_positive, iso_negative]

        #Do the actual rendering
        pov_obj.write(f'{label}.pov').render()

        os.remove(f'{label}.pov')
        os.remove(f'{label}.ini')

        if outfile != f'{label}.png': #move the output file to the desired location
            shutil.move(f'{label}.png', outfile)

    else: # use ASE renderer (low quality, does not draw bonds)
        write(outfile,
              atoms,
              format='png',
              radii = custom_settings.atomic_radius,
              rotation=rotations,
              colors=colors,
              maxwidth=width_res,
              scale=100)
