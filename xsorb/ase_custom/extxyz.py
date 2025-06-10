#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Enrico Pedretti
# Credit to original ASE code: https://wiki.fysik.dtu.dk/ase/

'''
Module to read and write extended XYZ files with custom labels
'''

import numpy as np
import ase.io.extxyz
from ase.constraints import FixAtoms, FixCartesian
from ase.io.extxyz import (XYZError, set_calc_and_arrays, key_val_str_to_dict, parse_properties,
    output_column_format, save_calc_results, voigt_6_to_full_3x3_stress)
from ase.utils import writer
from ase.io.espresso import label_to_symbol

from xsorb.ase_custom.atoms import AtomsCustom, extract_number_from_string


def _read_xyz_frame_custom(lines, natoms, properties_parser=key_val_str_to_dict,
                    nvec=0):
    """
    Custom version of ase.io.extxyz._read_xyz_frame to handle custom labels
    ---
    """

    # comment line
    line = next(lines).strip()
    if nvec > 0:
        info = {'comment': line}
    else:
        info = properties_parser(line) if line else {}

    pbc = None
    if 'pbc' in info:
        pbc = info.pop('pbc')
    elif 'Lattice' in info:
        # default pbc for extxyz file containing Lattice
        # is True in all directions
        pbc = [True, True, True]
    elif nvec > 0:
        # cell information given as pseudo-Atoms
        pbc = [False, False, False]

    cell = None
    if 'Lattice' in info:
        # NB: ASE cell is transpose of extended XYZ lattice
        cell = info['Lattice'].T
        del info['Lattice']
    elif nvec > 0:
        # cell information given as pseudo-Atoms
        cell = np.zeros((3, 3))

    if 'Properties' not in info:
        # Default set of properties is atomic symbols and positions only
        info['Properties'] = 'species:S:1:pos:R:3'
    properties, names, dtype, convs = parse_properties(info['Properties'])
    del info['Properties']

    data = []
    for _ in range(natoms):
        try:
            line = next(lines)
        except StopIteration:
            raise XYZError('ase.io.extxyz: Frame has {} atoms, expected {}'
                           .format(len(data), natoms))
        vals = line.split()
        row = tuple(conv(val) for conv, val in zip(convs, vals))
        data.append(row)

    try:
        data = np.array(data, dtype)
    except TypeError:
        raise XYZError('Badly formatted data '
                       'or end of file reached before end of frame')

    # Read VEC entries if present
    if nvec > 0:
        for ln in range(nvec):
            try:
                line = next(lines)
            except StopIteration:
                raise XYZError('ase.io.adfxyz: Frame has {} cell vectors, '
                               'expected {}'.format(len(cell), nvec))
            entry = line.split()

            if not entry[0].startswith('VEC'):
                raise XYZError(f'Expected cell vector, got {entry[0]}')

            try:
                n = int(entry[0][3:])
            except ValueError as e:
                raise XYZError('Expected VEC{}, got VEC{}'
                               .format(ln + 1, entry[0][3:])) from e
            if n != ln + 1:
                raise XYZError('Expected VEC{}, got VEC{}'
                               .format(ln + 1, n))

            cell[ln] = np.array([float(x) for x in entry[1:]])
            pbc[ln] = True
        if nvec != pbc.count(True):
            raise XYZError('Problem with number of cell vectors')
        pbc = tuple(pbc)

    arrays = {}
    for name in names:
        ase_name, cols = properties[name]
        if cols == 1:
            value = data[name]
        else:
            value = np.vstack([data[name + str(c)]
                               for c in range(cols)]).T
        arrays[ase_name] = value

    numbers = arrays.pop('numbers', None)
    symbols = arrays.pop('symbols', None)

    ##### CUSTOM PART #####
    if symbols is not None:
        symbols_plus_numbers = symbols.copy()
        symbols, tags = [], []
        for symbol_plus_number in symbols_plus_numbers:
            symbol = label_to_symbol(symbol_plus_number)
            number = extract_number_from_string(symbol_plus_number, symbol)
            symbols.append(symbol)
            tags.append(number)
        if not tags:
            tags = None
    else:
        tags=None

    atoms = AtomsCustom(numbers if numbers is not None else symbols,
                  positions=arrays.pop('positions', None),
                  charges=arrays.pop('initial_charges', None),
                  cell=cell,
                  pbc=pbc,
                  info=info,
                  tags=tags)
    ##### END CUSTOM PART #####

    # Read and set constraints
    if 'move_mask' in arrays:
        move_mask = arrays.pop('move_mask').astype(bool)
        if properties['move_mask'][1] == 3:
            cons = []
            for a in range(natoms):
                cons.append(FixCartesian(a, mask=~move_mask[a, :]))
            atoms.set_constraint(cons)
        elif properties['move_mask'][1] == 1:
            atoms.set_constraint(FixAtoms(mask=~move_mask))
        else:
            raise XYZError('Not implemented constraint')

    set_calc_and_arrays(atoms, arrays)
    return atoms


@writer
def write_xyz_custom(fileobj, images, comment='', columns=None,
              write_info=True,
              write_results=True, plain=False, vec_cell=False,
              custom_labels_as_symbols = True):
    """
    Custom version of ase.io.extxyz.write_xyz to handle custom_labels
    if custom_labels_as_symbols is True, custom_labels will be written as symbols
    ---


    Write output in extended XYZ format

    Optionally, specify which columns (arrays) to include in output,
    whether to write the contents of the `atoms.info` dict to the
    XYZ comment line (default is True), the results of any
    calculator attached to this Atoms. The `plain` argument
    can be used to write a simple XYZ file with no additional information.
    `vec_cell` can be used to write the cell vectors as additional
    pseudo-atoms. If `append` is set to True, the file is for append (mode `a`),
    otherwise it is overwritten (mode `w`).

    See documentation for :func:`read_xyz()` for further details of the extended
    XYZ file format.
    """
    if hasattr(images, 'get_positions'):
        images = [images]

    for atoms in images:
        natoms = len(atoms)

        if write_results:
            calculator = atoms.calc
            atoms = atoms.copy()

            save_calc_results(atoms, calculator, calc_prefix="")

            if atoms.info.get('stress', np.array([])).shape == (6,):
                atoms.info['stress'] = \
                    voigt_6_to_full_3x3_stress(atoms.info['stress'])

        if columns is None:
            fr_cols = (['symbols', 'positions']
                       + [key for key in atoms.arrays if
                          key not in ['symbols', 'positions', 'numbers',
                                      'species', 'pos']])
        else:
            fr_cols = columns[:]

        if vec_cell:
            plain = True

        if plain:
            fr_cols = ['symbols', 'positions']
            write_info = False
            write_results = False

        # Move symbols and positions to first two properties
        if 'symbols' in fr_cols:
            i = fr_cols.index('symbols')
            fr_cols[0], fr_cols[i] = fr_cols[i], fr_cols[0]

        if 'positions' in fr_cols:
            i = fr_cols.index('positions')
            fr_cols[1], fr_cols[i] = fr_cols[i], fr_cols[1]

        # Check first column "looks like" atomic symbols
        if fr_cols[0] in atoms.arrays:
            symbols = atoms.arrays[fr_cols[0]]
        else:
            symbols = [*atoms.symbols]

        if natoms > 0 and not isinstance(symbols[0], str):
            raise ValueError('First column must be symbols-like')

        # Check second column "looks like" atomic positions
        pos = atoms.arrays[fr_cols[1]]
        if pos.shape != (natoms, 3) or pos.dtype.kind != 'f':
            raise ValueError('Second column must be position-like')

        # if vec_cell add cell information as pseudo-atoms
        if vec_cell:
            nPBC = 0
            for i, b in enumerate(atoms.pbc):
                if not b:
                    continue
                nPBC += 1
                symbols.append('VEC' + str(nPBC))
                pos = np.vstack((pos, atoms.cell[i]))
            # add to natoms
            natoms += nPBC
            if pos.shape != (natoms, 3) or pos.dtype.kind != 'f':
                raise ValueError(
                    'Pseudo Atoms containing cell have bad coords')

        # Move mask
        if 'move_mask' in fr_cols:
            cnstr = images[0]._get_constraints()
            if len(cnstr) > 0:
                c0 = cnstr[0]
                if isinstance(c0, FixAtoms):
                    cnstr = np.ones((natoms,), dtype=bool)
                    for idx in c0.index:
                        cnstr[idx] = False  # cnstr: atoms that can be moved
                elif isinstance(c0, FixCartesian):
                    masks = np.ones((natoms, 3), dtype=bool)
                    for i in range(len(cnstr)):
                        idx = cnstr[i].index
                        masks[idx] = cnstr[i].mask
                    cnstr = ~masks  # cnstr: coordinates that can be moved
            else:
                fr_cols.remove('move_mask')

        # Collect data to be written out
        arrays = {}
        for column in fr_cols:
            if column == 'positions':
                arrays[column] = pos
            elif column in atoms.arrays:
                arrays[column] = atoms.arrays[column]
            elif column == 'symbols':
                arrays[column] = np.array(symbols)
            elif column == 'move_mask':
                arrays[column] = cnstr
            else:
                raise ValueError(f'Missing array "{column}"')

        ##### CUSTOM PART #####
        if custom_labels_as_symbols:
            arrays['symbols'] = atoms.custom_labels
        ##### END CUSTOM PART #####

        comm, ncols, dtype, fmt = output_column_format(atoms,
                                                       fr_cols,
                                                       arrays,
                                                       write_info)

        if plain or comment != '':
            # override key/value pairs with user-speficied comment string
            comm = comment.rstrip()
            if '\n' in comm:
                raise ValueError('Comment line should not have line breaks.')

        # Pack fr_cols into record array
        data = np.zeros(natoms, dtype)
        for column, ncol in zip(fr_cols, ncols):
            value = arrays[column]
            if ncol == 1:
                data[column] = np.squeeze(value)
            else:
                for c in range(ncol):
                    data[column + str(c)] = value[:, c]

        nat = natoms
        if vec_cell:
            nat -= nPBC
        # Write the output
        fileobj.write('%d\n' % nat)
        fileobj.write(f'{comm}\n')
        for i in range(natoms):
            fileobj.write(fmt % tuple(data[i]))


# Runtime patching
ase.io.extxyz._read_xyz_frame = _read_xyz_frame_custom
ase.io.extxyz.write_xyz = write_xyz_custom
