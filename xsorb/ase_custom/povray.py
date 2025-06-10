#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Enrico Pedretti
# Credit to original ASE code: https://wiki.fysik.dtu.dk/ase/

'''
Custom module to handle arrows in ase.io.pov, and to change fog and color styles
'''

import os
import sys

import numpy as np
import ase.io.utils
import ase.io.pov
from ase import Atoms
from ase.constraints import FixAtoms
from ase.data.colors import jmol_colors as default_colors
from ase.data import covalent_radii
from ase.io.pov import pa, pc


def POVRAYInit(self, cell, cell_vertices, positions, diameters, colors,
                 image_width, image_height, constraints=(), isosurfaces=[],
                 display=False, pause=True, transparent=True, canvas_width=None,
                 canvas_height=None, camera_dist=50., image_plane=None,
                 camera_type='orthographic', point_lights=[],
                 area_light=[(2., 3., 40.), 'White', .7, .7, 3, 3],
                 background='White', textures=None, transmittances=None,
                 depth_cueing=False, cue_density=5e-3, constant_fog_height=0.0,
                 celllinewidth=0.05, bondlinewidth=0.10, bondatoms=[],
                 exportconstraints=False,
                 arrows=None):
    """ Mod to ase.io.pov.POVRAY.__init__ with arrows """

    # attributes from initialization
    self.area_light = area_light
    self.background = background
    self.bondatoms = bondatoms
    self.bondlinewidth = bondlinewidth
    self.camera_dist = camera_dist
    self.camera_type = camera_type
    self.celllinewidth = celllinewidth
    self.cue_density = cue_density
    self.constant_fog_height = constant_fog_height
    self.depth_cueing = depth_cueing
    self.display = display
    self.exportconstraints = exportconstraints
    self.isosurfaces = isosurfaces
    self.pause = pause
    self.point_lights = point_lights
    self.textures = textures
    self.transmittances = transmittances
    self.transparent = transparent

    self.image_width = image_width
    self.image_height = image_height
    self.colors = colors
    self.cell = cell
    self.diameters = diameters
    self.arrows = arrows #### CUSTOM

    # calculations based on passed inputs

    #NOTE: fix, since PlottingVariables puts also some other strange coordinates besides those
    #of the atoms, so when calculating z0 the offset is wrong
    positions = positions[:len(diameters)] #### CUSTOM

    z0 = positions[:, 2].max()
    self.offset = (image_width / 2, image_height / 2, z0)
    self.positions = positions - self.offset

    if cell_vertices is not None:
        self.cell_vertices = cell_vertices - self.offset
        self.cell_vertices.shape = (2, 2, 2, 3)
    else:
        self.cell_vertices = None

    ratio = float(self.image_width) / self.image_height
    if canvas_width is None:
        if canvas_height is None:
            self.canvas_width = min(self.image_width * 15, 640)
            self.canvas_height = min(self.image_height * 15, 640)
        else:
            self.canvas_width = canvas_height * ratio
            self.canvas_height = canvas_height
    elif canvas_height is None:
        self.canvas_width = canvas_width
        self.canvas_height = self.canvas_width / ratio
    else:
        raise RuntimeError("Can't set *both* width and height!")

    # Distance to image plane from camera
    if image_plane is None:
        if self.camera_type == 'orthographic':
            self.image_plane = 1 - self.camera_dist
        else:
            self.image_plane = 0
    self.image_plane += self.camera_dist

    self.constrainatoms = []
    for c in constraints:
        if isinstance(c, FixAtoms):
            # self.constrainatoms.extend(c.index) # is this list-like?
            for n, i in enumerate(c.index):
                self.constrainatoms += [i]


@classmethod
def from_PlottingVariables(cls, pvars, **kwargs):
    '''
    Custom method to ase.io.pov.POVRAY.from_PlottingVariables
    to handle arrows

    '''
    cell = pvars.cell
    cell_vertices = pvars.cell_vertices
    if 'colors' in kwargs:
        colors = kwargs.pop('colors')
    else:
        colors = pvars.colors
    diameters = pvars.d
    image_height = pvars.h
    image_width = pvars.w
    positions = pvars.positions
    constraints = pvars.constraints
    arrows=pvars.arrows
    return cls(cell=cell, cell_vertices=cell_vertices, colors=colors, #pylint: disable=not-callable
            constraints=constraints, diameters=diameters,
            image_height=image_height, image_width=image_width,
            positions=positions, arrows=arrows, **kwargs) #### CUSTOM (arrows)


def PlottingVariablesInit(self, atoms, rotation='', show_unit_cell=2,
                 radii=None, bbox=None, colors=None, scale=20,
                 maxwidth=500, extra_offset=(0., 0.),
                 auto_bbox_size=1.05,
                 auto_image_plane_z='front_all',
                arrows_type=None):
    '''
    Custom init to ase.io.utils.PlottingVariables
    to handle arrows

    '''
    assert show_unit_cell in (0, 1, 2, 3)

    self.show_unit_cell = show_unit_cell
    self.numbers = atoms.get_atomic_numbers()
    self.maxwidth = maxwidth
    self.atoms = atoms
    # not used in PlottingVariables, keeping for legacy
    self.natoms = len(atoms)

    self.auto_bbox_size = auto_bbox_size
    self.auto_image_plane_z = auto_image_plane_z
    self.offset = np.zeros(3)
    self.extra_offset = np.array(extra_offset)

    self.constraints = atoms.constraints
    # extension for partial occupancies
    self.frac_occ = False
    self.tags = None
    self.occs = None

    if 'occupancy' in atoms.info:
        self.occs = atoms.info['occupancy']
        self.tags = atoms.get_tags()
        self.frac_occ = True

    # colors
    self.colors = colors
    if colors is None:
        ncolors = len(default_colors)
        self.colors = default_colors[self.numbers.clip(max=ncolors - 1)]

    # radius
    if radii is None:
        radii = covalent_radii[self.numbers]
    elif isinstance(radii, float):
        radii = covalent_radii[self.numbers] * radii
    else:
        radii = np.array(radii)

    self.radii = radii  # radius in Angstroms
    self.scale = scale  # Angstroms per cm

    self.set_rotation(rotation)
    self.update_image_plane_offset_and_size_from_structure(bbox=bbox)

    #### CUSTOM PART ####
    if arrows_type is not None:
        if arrows_type == 'forces':
            arrows = atoms.get_forces()
        elif arrows_type == 'magmoms':
            arrows = np.array([[0,0, magmom] for magmom in atoms.get_magnetic_moments()])
        else:
            raise ValueError(f'Unknown arrows type: {arrows_type}')
        arrows = np.dot(arrows, self.rotation)
        self.arrows = arrows
    else:
        self.arrows = None
    #### END CUSTOM PART ####


def write_ini(self, path):
    """
    Custom version of ase.io.pov.write_ini
    with Max_Image_Buffer_Memory=1024


    Write ini file."""

    ini_str = f"""\
Input_File_Name={path.with_suffix('.pov').name}
Output_to_File=True
Output_File_Type=N
Output_Alpha={'on' if self.transparent else 'off'}
; if you adjust Height, and width, you must preserve the ratio
; Width / Height = {self.canvas_width/self.canvas_height:f}
Width={self.canvas_width}
Height={self.canvas_height}
Antialias=True
Antialias_Threshold=0.1
Display={self.display}
Display_Gamma=2.2
Pause_When_Done={self.pause}
Verbose=False
Max_Image_Buffer_Memory=1024
"""
    with open(path, 'w', encoding=sys.getfilesystemencoding()) as fd:
        fd.write(ini_str)
    return path


def write_ini_old(self, path):
    """Write ini file."""

    ini_str = f"""\
Input_File_Name={path.with_suffix('.pov').name}
Output_to_File=True
Output_File_Type=N
Output_Alpha={'on' if self.transparent else 'off'}
; if you adjust Height, and width, you must preserve the ratio
; Width / Height = {self.canvas_width/self.canvas_height:f}
Width={self.canvas_width}
Height={self.canvas_height}
Antialias=True
Antialias_Threshold=0.1
Display={self.display}
Pause_When_Done={self.pause}
Verbose=False
Max_Image_Buffer_Memory=1024
"""
    with open(path, 'w', encoding=sys.getfilesystemencoding()) as fd:
        fd.write(ini_str)
    return path


def write_pov(self, path):
    """
    Custom version of ase.io.pov.write_pov with:
    - arrows
    - type 2 fog for depth cueing
    """

    point_lights = '\n'.join(f"light_source {{{pa(loc)} {pc(rgb)}}}"
                                for loc, rgb in self.point_lights)

    area_light = ''
    if self.area_light is not None:
        loc, color, width, height, nx, ny = self.area_light
        area_light += f"""\nlight_source {{{pa(loc)} {pc(color)}
area_light <{width:.2f}, 0, 0>, <0, {height:.2f}, 0>, {nx:n}, {ny:n}
adaptive 1 jitter}}"""

    fog = ''
    if self.depth_cueing and (self.cue_density >= 1e-4):
        # same way vmd does it
        if self.cue_density > 1e4:
            # larger does not make any sense
            dist = 1e-4
        else:
            dist = 2. / self.cue_density
        constant_fog_height = self.constant_fog_height if self.constant_fog_height is not None else 0.0
        fog += f'fog {{fog_type 2 distance {dist:.4f} up <0,0,1> fog_offset {constant_fog_height} fog_alt 0.1 '\
                f'color {pc(self.background)}}}'

    mat_style_keys = (f'#declare {k} = {v}'
                        for k, v in self.material_styles_dict.items())
    mat_style_keys = '\n'.join(mat_style_keys)

    # Draw unit cell
    cell_vertices = ''
    if self.cell_vertices is not None:
        for c in range(3):
            for j in ([0, 0], [1, 0], [1, 1], [0, 1]):
                p1 = self.cell_vertices[tuple(j[:c]) + (0,) + tuple(j[c:])]
                p2 = self.cell_vertices[tuple(j[:c]) + (1,) + tuple(j[c:])]

                distance = np.linalg.norm(p2 - p1)
                if distance < 1e-12:
                    continue

                cell_vertices += f'cylinder {{{pa(p1)}, {pa(p2)}, '\
                                    f'Rcell pigment {{Black}}}}\n'
                # all strings are f-strings for consistency
        cell_vertices = cell_vertices.strip('\n')

    # Draw atoms
    a = 0
    atoms = ''
    for loc, dia, col in zip(self.positions, self.diameters, self.colors):
        tex = 'ase3'
        trans = 0.
        if self.textures is not None:
            tex = self.textures[a]
        if self.transmittances is not None:
            trans = self.transmittances[a]
        atoms += f'atom({pa(loc)}, {dia/2.:.2f}, {pc(col)}, '\
                    f'{trans}, {tex}) // #{a:n}\n'
        a += 1
    atoms = atoms.strip('\n')

    # Draw atom bonds
    bondatoms = ''
    for pair in self.bondatoms:
        # Make sure that each pair has 4 componets: a, b, offset,
        #                                           bond_order, bond_offset
        # a, b: atom index to draw bond
        # offset: original meaning to make offset for mid-point.
        # bond_oder: if not supplied, set it to 1 (single bond).
        #            It can be  1, 2, 3, corresponding to single,
        #            double, triple bond
        # bond_offset: displacement from original bond position.
        #              Default is (bondlinewidth, bondlinewidth, 0)
        #              for bond_order > 1.
        if len(pair) == 2:
            a, b = pair
            offset = (0, 0, 0)
            bond_order = 1
            bond_offset = (0, 0, 0)
        elif len(pair) == 3:
            a, b, offset = pair
            bond_order = 1
            bond_offset = (0, 0, 0)
        elif len(pair) == 4:
            a, b, offset, bond_order = pair
            bond_offset = (self.bondlinewidth, self.bondlinewidth, 0)
        elif len(pair) > 4:
            a, b, offset, bond_order, bond_offset = pair
        else:
            raise RuntimeError('Each list in bondatom must have at least '
                                '2 entries. Error at %s' % pair)

        if len(offset) != 3:
            raise ValueError('offset must have 3 elements. '
                                'Error at %s' % pair)
        if len(bond_offset) != 3:
            raise ValueError('bond_offset must have 3 elements. '
                                'Error at %s' % pair)
        if bond_order not in [0, 1, 2, 3]:
            raise ValueError('bond_order must be either 0, 1, 2, or 3. '
                                'Error at %s' % pair)

        # Up to here, we should have all a, b, offset, bond_order,
        # bond_offset for all bonds.

        # Rotate bond_offset so that its direction is 90 deg. off the bond
        # Utilize Atoms object to rotate
        if bond_order > 1 and np.linalg.norm(bond_offset) > 1.e-9:
            tmp_atoms = Atoms('H3')
            tmp_atoms.set_cell(self.cell)
            tmp_atoms.set_positions([
                self.positions[a],
                self.positions[b],
                self.positions[b] + np.array(bond_offset),
            ])
            tmp_atoms.center()
            tmp_atoms.set_angle(0, 1, 2, 90)
            bond_offset = tmp_atoms[2].position - tmp_atoms[1].position

        R = np.dot(offset, self.cell)
        mida = 0.5 * (self.positions[a] + self.positions[b] + R)
        midb = 0.5 * (self.positions[a] + self.positions[b] - R)
        if self.textures is not None:
            texa = self.textures[a]
            texb = self.textures[b]
        else:
            texa = texb = 'ase3'

        if self.transmittances is not None:
            transa = self.transmittances[a]
            transb = self.transmittances[b]
        else:
            transa = transb = 0.

        # draw bond, according to its bond_order.
        # bond_order == 0: No bond is plotted
        # bond_order == 1: use original code
        # bond_order == 2: draw two bonds, one is shifted by bond_offset/2,
        #                  and another is shifted by -bond_offset/2.
        # bond_order == 3: draw two bonds, one is shifted by bond_offset,
        #                  and one is shifted by -bond_offset, and the
        #                  other has no shift.
        # To shift the bond, add the shift to the first two coordinate in
        # write statement.

        posa = self.positions[a]
        posb = self.positions[b]
        cola = self.colors[a]
        colb = self.colors[b]

        if bond_order == 1:
            draw_tuples = (
                (posa, mida, cola, transa, texa),
                (posb, midb, colb, transb, texb))

        elif bond_order == 2:
            bs = [x / 2 for x in bond_offset]
            draw_tuples = (
                (posa - bs, mida - bs, cola, transa, texa),
                (posb - bs, midb - bs, colb, transb, texb),
                (posa + bs, mida + bs, cola, transa, texa),
                (posb + bs, midb + bs, colb, transb, texb))

        elif bond_order == 3:
            bs = bond_offset
            draw_tuples = (
                (posa, mida, cola, transa, texa),
                (posb, midb, colb, transb, texb),
                (posa + bs, mida + bs, cola, transa, texa),
                (posb + bs, midb + bs, colb, transb, texb),
                (posa - bs, mida - bs, cola, transa, texa),
                (posb - bs, midb - bs, colb, transb, texb))

        bondatoms += ''.join(f'cylinder {{{pa(p)}, '
                                f'{pa(m)}, Rbond texture{{pigment '
                                f'{{color {pc(c)} '
                                f'transmit {tr}}} finish{{{tx}}}}}}}\n'
                                for p, m, c, tr, tx in
                                draw_tuples)

    bondatoms = bondatoms.strip('\n')

    # Draw constraints if requested
    constraints = ''
    if self.exportconstraints:
        for a in self.constrainatoms:
            dia = self.diameters[a]
            loc = self.positions[a]
            trans = 0.0
            if self.transmittances is not None:
                trans = self.transmittances[a]
            constraints += f'constrain({pa(loc)}, {dia / 2.:.2f}, Black, '\
                f'{trans}, {tex}) // #{a:n} \n'
    constraints = constraints.strip('\n')

    #### BEGIN CUSTOM: handle arrows
    # Draw arrows
    arrows = ''
    if self.arrows is not None:
        maxlength = max([np.linalg.norm(arrow) for arrow in self.arrows])
        for pos, arrow, diam in zip(self.positions, self.arrows, self.diameters):
            modulus = np.linalg.norm(arrow)
            normalized_arrow = arrow / maxlength
            if modulus/maxlength > 0.1: # avoid degenerate primitives
                cylinder_pos_dw = pos - 0.8*normalized_arrow
                cylinder_pos_up = pos + 0.7*normalized_arrow
                arrows += f'cylinder {{{pa(cylinder_pos_dw)}, '+\
                                        f'{pa(cylinder_pos_up)}, 0.1 texture{{pigment '+\
                                        f'{{color {pc([1,0,0])} '+\
                                        f'transmit 0.0}} finish{{ase3}}}}}}\n'
                cone_pos = pos + 0.7*normalized_arrow
                arrows += f'cone {{{pa(cone_pos)}, 0.2'+\
                                        f'{pa(cone_pos + 0.3*arrow/modulus)}, 0.0 texture{{pigment '+\
                                        f'{{color {pc([1,0,0])} '+\
                                        f'transmit 0.0}} finish{{ase3}}}}}}\n'
    #### END CUSTOM

    pov = f"""#version 3.6;
#include "colors.inc"
#include "finish.inc"

global_settings {{assumed_gamma 2.2 max_trace_level 6}}
background {{{pc(self.background)}{' transmit 1.0' if self.transparent else ''}}}
camera {{{self.camera_type}
right -{self.image_width:.2f}*x up {self.image_height:.2f}*y
direction {self.image_plane:.2f}*z
location <0,0,{self.camera_dist:.2f}> look_at <0,0,0>}}
{point_lights}
{area_light if area_light != '' else '// no area light'}
{fog if fog != '' else '// no fog'}
{mat_style_keys}
#declare Rcell = {self.celllinewidth:.3f};
#declare Rbond = {self.bondlinewidth:.3f};

#macro atom(LOC, R, COL, TRANS, FIN)
sphere{{LOC, R texture{{pigment{{color COL transmit TRANS}} finish{{FIN}}}}}}
#end
#macro constrain(LOC, R, COL, TRANS FIN)
union{{torus{{R, Rcell rotate 45*z texture{{pigment{{color COL transmit TRANS}} finish{{FIN}}}}}}
    torus{{R, Rcell rotate -45*z texture{{pigment{{color COL transmit TRANS}} finish{{FIN}}}}}}
    translate LOC}}
#end

{cell_vertices if cell_vertices != '' else '// no cell vertices'}
{atoms}
{bondatoms}
{constraints if constraints != '' else '// no constraints'}
{arrows if arrows != '' else '// no arrows'}
"""  # noqa: E501

    with open(path, 'w') as fd:
        fd.write(pov)

    return path



# Runtime patching
if 'POVRAY_OLD_STYLE' in os.environ:
    ase.io.pov.POVRAY.write_ini = write_ini_old
    ase.io.pov.POVRAY.material_styles_dict = ase.io.pov.POVRAY.material_styles_dict_old

else:
    ase.io.pov.POVRAY.write_ini = write_ini
    ase.io.pov.POVRAY.material_styles_dict['pale']=('finish {ambient 0.9 diffuse 0.30 roughness 0.001}')

ase.io.pov.POVRAY.__init__ = POVRAYInit
ase.io.pov.POVRAY.write_pov = write_pov
ase.io.pov.POVRAY.from_PlottingVariables = from_PlottingVariables
ase.io.utils.PlottingVariables.__init__ = PlottingVariablesInit
