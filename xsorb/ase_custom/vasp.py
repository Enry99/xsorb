#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Enrico Pedretti
# Credit to original ASE code: https://wiki.fysik.dtu.dk/ase/


'''
Mod to ase.io.vasp.read_vasp to add resorting for poscar if ase-sort.dat is preset
'''

from pathlib import Path
import re
import sys

import numpy as np
import ase.io.vasp
from ase import Atoms
from ase.utils import reader
from ase.io.vasp import (get_atomtypes_from_formula, atomtypes_outpot,
    set_constraints, parse_poscar_scaling_factor)

@reader
def read_vasp_configuration_custom(fd):
    """
    Custom mod to ase.io.vasp.read_vasp_configuration to add resorting
    for POSCAR if ase-sort.dat is present.
    ----

    Read common POSCAR/CONTCAR/CHGCAR/CHG quantities and return Atoms.
    """
    from ase.data import chemical_symbols

    # The first line is in principle a comment line, however in VASP
    # 4.x a common convention is to have it contain the atom symbols,
    # eg. "Ag Ge" in the same order as later in the file (and POTCAR
    # for the full vasp run). In the VASP 5.x format this information
    # is found on the fifth line. Thus we save the first line and use
    # it in case we later detect that we're reading a VASP 4.x format
    # file.
    line1 = fd.readline()

    scale = parse_poscar_scaling_factor(fd.readline())

    # Now the lattice vectors
    cell = np.array([fd.readline().split()[:3] for _ in range(3)], dtype=float)
    # Negative scaling factor corresponds to the cell volume.
    if scale[0] < 0.0:
        scale = np.cbrt(-1.0 * scale / np.linalg.det(cell))
    cell *= scale  # This works for both one and three scaling factors.

    # Number of atoms. Again this must be in the same order as
    # in the first line
    # or in the POTCAR or OUTCAR file
    atom_symbols = []
    numofatoms = fd.readline().split()
    # Check whether we have a VASP 4.x or 5.x format file. If the
    # format is 5.x, use the fifth line to provide information about
    # the atomic symbols.
    vasp5 = False
    try:
        int(numofatoms[0])
    except ValueError:
        vasp5 = True
        atomtypes = numofatoms
        numofatoms = fd.readline().split()

    # check for comments in numofatoms line and get rid of them if necessary
    commentcheck = np.array(['!' in s for s in numofatoms])
    if commentcheck.any():
        # only keep the elements up to the first including a '!':
        numofatoms = numofatoms[:np.arange(len(numofatoms))[commentcheck][0]]

    if not vasp5:
        # Split the comment line (first in the file) into words and
        # try to compose a list of chemical symbols
        from ase.formula import Formula
        atomtypes = []
        for word in line1.split():
            word_without_delims = re.sub(r"-|_|,|\.|=|[0-9]|^", "", word)
            if len(word_without_delims) < 1:
                continue
            try:
                atomtypes.extend(list(Formula(word_without_delims)))
            except ValueError:
                # print(atomtype, e, 'is comment')
                pass
        # Now the list of chemical symbols atomtypes must be formed.
        # For example: atomtypes = ['Pd', 'C', 'O']

        numsyms = len(numofatoms)
        if len(atomtypes) < numsyms:
            # First line in POSCAR/CONTCAR didn't contain enough symbols.

            # Sometimes the first line in POSCAR/CONTCAR is of the form
            # "CoP3_In-3.pos". Check for this case and extract atom types
            if len(atomtypes) == 1 and '_' in atomtypes[0]:
                atomtypes = get_atomtypes_from_formula(atomtypes[0])
            else:
                atomtypes = atomtypes_outpot(fd.name, numsyms)
        else:
            try:
                for atype in atomtypes[:numsyms]:
                    if atype not in chemical_symbols:
                        raise KeyError
            except KeyError:
                atomtypes = atomtypes_outpot(fd.name, numsyms)

    for i, num in enumerate(numofatoms):
        numofatoms[i] = int(num)
        atom_symbols.extend(numofatoms[i] * [atomtypes[i]])

    # Check if Selective dynamics is switched on
    sdyn = fd.readline()
    selective_dynamics = sdyn[0].lower() == 's'

    # Check if atom coordinates are cartesian or direct
    if selective_dynamics:
        ac_type = fd.readline()
    else:
        ac_type = sdyn
    cartesian = ac_type[0].lower() in ['c', 'k']
    tot_natoms = sum(numofatoms)
    atoms_pos = np.empty((tot_natoms, 3))
    if selective_dynamics:
        selective_flags = np.empty((tot_natoms, 3), dtype=bool)
    for atom in range(tot_natoms):
        ac = fd.readline().split()
        atoms_pos[atom] = [float(_) for _ in ac[0:3]]
        if selective_dynamics:
            selective_flags[atom] = [_ == 'F' for _ in ac[3:6]]

    atoms = Atoms(symbols=atom_symbols, cell=cell, pbc=True)
    if cartesian:
        atoms_pos *= scale
        atoms.set_positions(atoms_pos)
    else:
        atoms.set_scaled_positions(atoms_pos)
    if selective_dynamics:
        set_constraints(atoms, selective_flags)

    #### CUSTOM PART ####
    # Resort if ase-sort.dat is present
    sortfile = f'{Path(fd.name).parent}/ase-sort.dat'
    if Path(sortfile).exists():
        resort_list = []
        with open(sortfile, 'r') as ff:
            for line in ff:
                sort, resort = line.split()
                resort_list.append(int(resort))

        atoms = atoms[resort_list]
    #### END CUSTOM PART ####

    return atoms

# Runtime patching
ase.io.vasp.read_vasp_configuration = read_vasp_configuration_custom
