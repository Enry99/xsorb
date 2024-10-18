#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Enrico Pedretti
# Credit to original ASE code: https://wiki.fysik.dtu.dk/ase/

'''
Module to customize the ASE interface to Quantum Espresso.
It includes custom labels for atoms, and the ability to pass through pseudopotentials.
'''

import warnings
import re

import numpy as np
import ase.io.espresso
from ase.io.espresso import read_fortran_namelist, get_cell_parameters, ibrav_error_message, \
    get_atomic_positions, get_atomic_species, convert_constraint_flags, label_to_symbol, units, \
    kspacing_to_grid, kpts2sizeandoffsets, kpts2ndarray, Namelist, parse_position_line,\
    _PW_START,_PW_END,_PW_CELL,_PW_POS,_PW_MAGMOM,_PW_FORCE, _PW_TOTEN,_PW_STRESS,_PW_FERMI,\
    _PW_HIGHEST_OCCUPIED,_PW_HIGHEST_OCCUPIED_LOWEST_FREE,_PW_KPTS,_PW_BANDS,_PW_BANDSTRUCTURE,\
    _PW_DIPOLE, _PW_DIPOLE_DIRECTION, kpoint_convert
from ase.calculators.singlepoint import SinglePointDFTCalculator,SinglePointKPoint
from ase.utils import reader
from ase import Atom
from ase.constraints import FixAtoms, FixCartesian

from xsorb.ase_custom.atoms import AtomsCustom


@reader
def read_espresso_in_custom(fileobj):
    """Patched version of ase.io.espresso.read_espresso_in to include custom labels,
    e.g. 'Fe1', 'H2'.
    """
    # parse namelist section and extract remaining lines
    data, card_lines = read_fortran_namelist(fileobj)

    # get the cell if ibrav=0
    if 'system' not in data:
        raise KeyError('Required section &SYSTEM not found.')
    elif 'ibrav' not in data['system']:
        raise KeyError('ibrav is required in &SYSTEM')
    elif data['system']['ibrav'] == 0:
        # celldm(1) is in Bohr, A is in angstrom. celldm(1) will be
        # used even if A is also specified.
        if 'celldm(1)' in data['system']:
            alat = data['system']['celldm(1)'] * units['Bohr']
        elif 'A' in data['system']:
            alat = data['system']['A']
        else:
            alat = None
        cell, _ = get_cell_parameters(card_lines, alat=alat)
    else:
        raise ValueError(ibrav_error_message)

    # species_info holds some info for each element
    species_card = get_atomic_species(
        card_lines, n_species=data['system']['ntyp'])
    species_info = {}
    for ispec, (label, weight, pseudo) in enumerate(species_card):
        symbol = label   # NOTE: changed to handle custom labels

        # starting_magnetization is in fractions of valence electrons
        magnet_key = f"starting_magnetization({ispec + 1})"
        magmom = data["system"].get(magnet_key, 0.0)
        species_info[symbol] = {"weight": weight, "pseudo": pseudo,
                                "magmom": magmom}


    positions_card = get_atomic_positions(
        card_lines, n_atoms=data['system']['nat'], cell=cell, alat=alat)

    symbols = [label_to_symbol(position[0]) for position in positions_card]
    custom_labels = [position[0] for position in positions_card]
    positions = [position[1] for position in positions_card]
    constraint_flags = [position[2] for position in positions_card]
    magmoms = [species_info[symbol]["magmom"] for symbol in symbols]


    # TODO: put more info into the atoms object
    # e.g magmom, forces.
    atoms = AtomsCustom(symbols=symbols, positions=positions, cell=cell, pbc=True,
                  magmoms=magmoms, custom_labels=custom_labels)
    atoms.set_constraint(convert_constraint_flags(constraint_flags))

    return atoms


def format_atom_position(atom, crystal_coordinates, custom_label, mask='', tidx=None):
    """
    Custom version to handle custom labels

    Format one line of atomic positions in
    Quantum ESPRESSO ATOMIC_POSITIONS card.

    >>> for atom in make_supercell(bulk('Li', 'bcc'), np.ones(3)-np.eye(3)):
    >>>     format_atom_position(atom, True)
    Li 0.0000000000 0.0000000000 0.0000000000
    Li 0.5000000000 0.5000000000 0.5000000000

    Parameters
    ----------
    atom : Atom
        A structure that has symbol and [position | (a, b, c)].
    crystal_coordinates: bool
        Whether the atomic positions should be written to the QE input file in
        absolute (False, default) or relative (crystal) coordinates (True).
    mask, optional : str
        String of ndim=3 0 or 1 for constraining atomic positions.
    tidx, optional : int
        Magnetic type index.

    Returns
    -------
    atom_line : str
        Input line for atom position
    """
    if crystal_coordinates:
        coords = [atom.a, atom.b, atom.c]
    else:
        coords = atom.position
    line_fmt = f'{custom_label}'
    inps = dict(atom=atom)
    if tidx is not None:
        line_fmt += '{tidx}'
        inps["tidx"] = tidx
    line_fmt += ' {coords[0]:.10f} {coords[1]:.10f} {coords[2]:.10f} '
    inps["coords"] = coords
    line_fmt += ' ' + mask + '\n'
    astr = line_fmt.format(**inps)
    return astr


def write_espresso_in_custom(fd, atoms, input_data=None, pseudopotentials=None,
                      kspacing=None, kpts=None, koffset=(0, 0, 0),
                      crystal_coordinates=False, additional_cards=None,
                      **kwargs):
    """
    Patched version of ase.io.espresso.write_espresso_in,
    which simply passes through the pseudopotentials, leaving the order unchanged,
    and keeps the labels of the atoms
    """


    # Convert to a namelist to make working with parameters much easier
    # Note that the name ``input_data`` is chosen to prevent clash with
    # ``parameters`` in Calculator objects
    input_parameters = Namelist(input_data)
    input_parameters.to_nested('pw', **kwargs)

    # Convert ase constraints to QE constraints
    # Nx3 array of force multipliers matches what QE uses
    # Do this early so it is available when constructing the atoms card
    moved = np.ones((len(atoms), 3), dtype=bool)
    for constraint in atoms.constraints:
        if isinstance(constraint, FixAtoms):
            moved[constraint.index] = False
        elif isinstance(constraint, FixCartesian):
            moved[constraint.index] = ~constraint.mask
        else:
            warnings.warn(f'Ignored unknown constraint {constraint}')
    masks = []
    for atom in atoms:
        # only inclued mask if something is fixed
        if not all(moved[atom.index]):
            mask = ' {:d} {:d} {:d}'.format(*moved[atom.index])
        else:
            mask = ''
        masks.append(mask)

    # Species info holds the information on the pseudopotential and
    # associated for each element
    if pseudopotentials is None:
        pseudopotentials = {}
    species_info = {}
    for species in set(atoms.get_chemical_symbols()):
        # Look in all possible locations for the pseudos and try to figure
        # out the number of valence electrons
        pseudo = pseudopotentials[species]
        species_info[species] = {'pseudo': pseudo}

    # Convert atoms into species.
    # Each different magnetic moment needs to be a separate type even with
    # the same pseudopotential (e.g. an up and a down for AFM).
    # if any magmom are > 0 or nspin == 2 then use species labels.
    # Rememeber: magnetisation uses 1 based indexes
    atomic_species = {}
    atomic_species_str = []
    atomic_positions_str = []

    #do not force spin on, simply passthrough
    # nspin = input_parameters['system'].get('nspin', 1)  # 1 is the default
    # noncolin = input_parameters['system'].get('noncolin', False)
    # rescale_magmom_fac = kwargs.get('rescale_magmom_fac', 1.0)
    # if any(atoms.get_initial_magnetic_moments()):
    #     if nspin == 1 and not noncolin:
    #         # Force spin on
    #         input_parameters['system']['nspin'] = 2
    #         nspin = 2

    # if nspin == 2 or noncolin:
    #     # Magnetic calculation on
    #     for atom, mask, magmom in zip(
    #             atoms, masks, atoms.get_initial_magnetic_moments()):
    #         if (atom.symbol, magmom) not in atomic_species:
    #             # for qe version 7.2 or older magmon must be rescale by
    #             # about a factor 10 to assume sensible values
    #             # since qe-v7.3 magmom values will be provided unscaled
    #             fspin = float(magmom) / rescale_magmom_fac
    #             # Index in the atomic species list
    #             sidx = len(atomic_species) + 1
    #             # Index for that atom type; no index for first one
    #             tidx = sum(atom.symbol == x[0] for x in atomic_species) or ' '
    #             atomic_species[(atom.symbol, magmom)] = (sidx, tidx)
    #             # Add magnetization to the input file
    #             mag_str = f"starting_magnetization({sidx})"
    #             input_parameters['system'][mag_str] = fspin
    #             species_pseudo = species_info[atom.symbol]['pseudo']
    #             atomic_species_str.append(
    #                 f"{atom.symbol}{tidx} {atom.mass} {species_pseudo}\n")
    #         # lookup tidx to append to name
    #         sidx, tidx = atomic_species[(atom.symbol, magmom)]
    #         # construct line for atomic positions
    #         atomic_positions_str.append(
    #             format_atom_position(
    #                 atom, crystal_coordinates, mask=mask, tidx=tidx)
    #         )
    # else:
    # Do nothing about magnetisation

    #CHANGED THIS PART FOR PASSTHROUGH AND HANDLE CUSTOM LABELS
    for label, pseudo in pseudopotentials.items():
        atomic_species[label] = True  # just a placeholder
        atomic_species_str.append(
                '{label} {mass} {pseudo}\n'.format(
                    label=label,
                    mass=Atom(label_to_symbol(label)).mass,
                    pseudo=pseudo))

    for atom, mask, custom_label in zip(atoms, masks, atoms.get_custom_labels()):
        # construct line for atomic positions
        atomic_positions_str.append(
            format_atom_position(atom, crystal_coordinates, custom_label, mask=mask)
        )
    ##########################################################

    # Add computed parameters
    # different magnetisms means different types
    input_parameters['system']['ntyp'] = len(atomic_species)
    input_parameters['system']['nat'] = len(atoms)

    # Use cell as given or fit to a specific ibrav
    if 'ibrav' in input_parameters['system']:
        ibrav = input_parameters['system']['ibrav']
        if ibrav != 0:
            raise ValueError(ibrav_error_message)
    else:
        # Just use standard cell block
        input_parameters['system']['ibrav'] = 0

    # Construct input file into this
    pwi = input_parameters.to_string(list_form=True)

    # Pseudopotentials
    pwi.append('ATOMIC_SPECIES\n')
    pwi.extend(atomic_species_str)
    pwi.append('\n')

    # KPOINTS - add a MP grid as required
    if kspacing is not None:
        kgrid = kspacing_to_grid(atoms, kspacing)
    elif kpts is not None:
        if isinstance(kpts, dict) and 'path' not in kpts:
            kgrid, shift = kpts2sizeandoffsets(atoms=atoms, **kpts)
            koffset = []
            for i, x in enumerate(shift):
                assert x == 0 or abs(x * kgrid[i] - 0.5) < 1e-14
                koffset.append(0 if x == 0 else 1)
        else:
            kgrid = kpts
    else:
        kgrid = "gamma"

    # True and False work here and will get converted by ':d' format
    if isinstance(koffset, int):
        koffset = (koffset, ) * 3

    # BandPath object or bandpath-as-dictionary:
    if isinstance(kgrid, dict) or hasattr(kgrid, 'kpts'):
        pwi.append('K_POINTS crystal_b\n')
        assert hasattr(kgrid, 'path') or 'path' in kgrid
        kgrid = kpts2ndarray(kgrid, atoms=atoms)
        pwi.append(f'{len(kgrid)}\n')
        for k in kgrid:
            pwi.append(f"{k[0]:.14f} {k[1]:.14f} {k[2]:.14f} 0\n")
        pwi.append('\n')
    elif isinstance(kgrid, str) and (kgrid == "gamma"):
        pwi.append('K_POINTS gamma\n')
        pwi.append('\n')
    else:
        pwi.append('K_POINTS automatic\n')
        pwi.append(f"{kgrid[0]} {kgrid[1]} {kgrid[2]} "
                   f" {koffset[0]:d} {koffset[1]:d} {koffset[2]:d}\n")
        pwi.append('\n')

    # CELL block, if required
    if input_parameters['SYSTEM']['ibrav'] == 0:
        pwi.append('CELL_PARAMETERS angstrom\n')
        pwi.append('{cell[0][0]:.14f} {cell[0][1]:.14f} {cell[0][2]:.14f}\n'
                   '{cell[1][0]:.14f} {cell[1][1]:.14f} {cell[1][2]:.14f}\n'
                   '{cell[2][0]:.14f} {cell[2][1]:.14f} {cell[2][2]:.14f}\n'
                   ''.format(cell=atoms.cell))
        pwi.append('\n')

    # Positions - already constructed, but must appear after namelist
    if crystal_coordinates:
        pwi.append('ATOMIC_POSITIONS crystal\n')
    else:
        pwi.append('ATOMIC_POSITIONS angstrom\n')
    pwi.extend(atomic_positions_str)
    pwi.append('\n')

    # DONE!
    fd.write(''.join(pwi))

    if additional_cards:
        if isinstance(additional_cards, list):
            additional_cards = "\n".join(additional_cards)
            additional_cards += "\n"

        fd.write(additional_cards)


def parse_pwo_start_custom(lines, index=0):
    """
    Custom version of ase.io.espresso.parse_pwo_start
    that handles custom labels

    Parse Quantum ESPRESSO calculation info from lines,
    starting from index. Return a dictionary containing extracted
    information.

    - `celldm(1)`: lattice parameters (alat)
    - `cell`: unit cell in Angstrom
    - `symbols`: element symbols for the structure
    - `positions`: cartesian coordinates of atoms in Angstrom
    - `atoms`: an `ase.Atoms` object constructed from the extracted data

    Parameters
    ----------
    lines : list[str]
        Contents of PWSCF output file.
    index : int
        Line number to begin parsing. Only first calculation will
        be read.

    Returns
    -------
    info : dict
        Dictionary of calculation parameters, including `celldm(1)`, `cell`,
        `symbols`, `positions`, `atoms`.

    Raises
    ------
    KeyError
        If interdependent values cannot be found (especially celldm(1))
        an error will be raised as other quantities cannot then be
        calculated (e.g. cell and positions).
    """

    info = {}

    custom_labels = []

    for idx, line in enumerate(lines[index:], start=index):
        if 'celldm(1)' in line:
            # celldm(1) has more digits than alat!!
            info['celldm(1)'] = float(line.split()[1]) * units['Bohr']
            info['alat'] = info['celldm(1)']
        elif 'number of atoms/cell' in line:
            info['nat'] = int(line.split()[-1])
        elif 'number of atomic types' in line:
            info['ntyp'] = int(line.split()[-1])
        elif 'crystal axes:' in line:
            info['cell'] = info['celldm(1)'] * np.array([
                [float(x) for x in lines[idx + 1].split()[3:6]],
                [float(x) for x in lines[idx + 2].split()[3:6]],
                [float(x) for x in lines[idx + 3].split()[3:6]]])
        elif 'positions (alat units)' in line:
            info['symbols'], info['positions'] = [], []

            for at_line in lines[idx + 1:idx + 1 + info['nat']]:
                sym, x, y, z = parse_position_line(at_line)
                info['symbols'].append(label_to_symbol(sym))
                custom_labels.append(sym)
                info['positions'].append([x * info['celldm(1)'],
                                          y * info['celldm(1)'],
                                          z * info['celldm(1)']])
            # This should be the end of interesting info.
            # Break here to avoid dealing with large lists of kpoints.
            # Will need to be extended for DFTCalculator info.
            break

    # Make atoms for convenience
    info['atoms'] = AtomsCustom(symbols=info['symbols'],
                          positions=info['positions'],
                          cell=info['cell'], pbc=True, custom_labels=custom_labels)

    return info

#TODO: in the future, fill in with the one from MDAnalysis,
# which hanldes the restarts.
@reader
def read_espresso_out_custom(fileobj, index=-1, results_required=True):
    """
    Custom version of ase.io.espresso.read_espresso_out
    that handles custom labels


    Reads Quantum ESPRESSO output files.

    The atomistic configurations as well as results (energy, force, stress,
    magnetic moments) of the calculation are read for all configurations
    within the output file.

    Will probably raise errors for broken or incomplete files.

    Parameters
    ----------
    fileobj : file|str
        A file like object or filename
    index : slice
        The index of configurations to extract.
    results_required : bool
        If True, atomistic configurations that do not have any
        associated results will not be included. This prevents double
        printed configurations and incomplete calculations from being
        returned as the final configuration with no results data.

    Yields
    ------
    structure : Atoms
        The next structure from the index slice. The Atoms has a
        SinglePointCalculator attached with any results parsed from
        the file.


    """
    # work with a copy in memory for faster random access
    pwo_lines = fileobj.readlines()

    # TODO: index -1 special case?
    # Index all the interesting points
    indexes = {
        _PW_START: [],
        _PW_END: [],
        _PW_CELL: [],
        _PW_POS: [],
        _PW_MAGMOM: [],
        _PW_FORCE: [],
        _PW_TOTEN: [],
        _PW_STRESS: [],
        _PW_FERMI: [],
        _PW_HIGHEST_OCCUPIED: [],
        _PW_HIGHEST_OCCUPIED_LOWEST_FREE: [],
        _PW_KPTS: [],
        _PW_BANDS: [],
        _PW_BANDSTRUCTURE: [],
    }

    for idx, line in enumerate(pwo_lines):
        for identifier in indexes:
            if identifier in line:
                indexes[identifier].append(idx)

    # Configurations are either at the start, or defined in ATOMIC_POSITIONS
    # in a subsequent step. Can deal with concatenated output files.
    all_config_indexes = sorted(indexes[_PW_START] +
                                indexes[_PW_POS])

    # Slice only requested indexes
    # setting results_required argument stops configuration-only
    # structures from being returned. This ensures the [-1] structure
    # is one that has results. Two cases:
    # - SCF of last configuration is not converged, job terminated
    #   abnormally.
    # - 'relax' and 'vc-relax' re-prints the final configuration but
    #   only 'vc-relax' recalculates.
    if results_required:
        results_indexes = sorted(indexes[_PW_TOTEN] + indexes[_PW_FORCE] +
                                 indexes[_PW_STRESS] + indexes[_PW_MAGMOM] +
                                 indexes[_PW_BANDS] +
                                 indexes[_PW_BANDSTRUCTURE])

        # Prune to only configurations with results data before the next
        # configuration
        results_config_indexes = []
        for config_index, config_index_next in zip(
                all_config_indexes,
                all_config_indexes[1:] + [len(pwo_lines)]):
            if any(config_index < results_index < config_index_next
                    for results_index in results_indexes):
                results_config_indexes.append(config_index)

        # slice from the subset
        image_indexes = results_config_indexes[index]
    else:
        image_indexes = all_config_indexes[index]

    # Extract initialisation information each time PWSCF starts
    # to add to subsequent configurations. Use None so slices know
    # when to fill in the blanks.
    pwscf_start_info = {idx: None for idx in indexes[_PW_START]}

    for image_index in image_indexes:
        # Find the nearest calculation start to parse info. Needed in,
        # for example, relaxation where cell is only printed at the
        # start.
        if image_index in indexes[_PW_START]:
            prev_start_index = image_index
        else:
            # The greatest start index before this structure
            prev_start_index = [idx for idx in indexes[_PW_START]
                                if idx < image_index][-1]

        # add structure to reference if not there
        if pwscf_start_info[prev_start_index] is None:
            pwscf_start_info[prev_start_index] = parse_pwo_start_custom(
                pwo_lines, prev_start_index)

        # Get the bounds for information for this structure. Any associated
        # values will be between the image_index and the following one,
        # EXCEPT for cell, which will be 4 lines before if it exists.
        for next_index in all_config_indexes:
            if next_index > image_index:
                break
        else:
            # right to the end of the file
            next_index = len(pwo_lines)

        # Get the structure
        # Use this for any missing data
        prev_structure = pwscf_start_info[prev_start_index]['atoms']
        if image_index in indexes[_PW_START]:
            structure = prev_structure.copy()  # parsed from start info
        else:
            if _PW_CELL in pwo_lines[image_index - 5]:
                # CELL_PARAMETERS would be just before positions if present
                cell, cell_alat = get_cell_parameters(
                    pwo_lines[image_index - 5:image_index])
            else:
                cell = prev_structure.cell
                cell_alat = pwscf_start_info[prev_start_index]['alat']

            # give at least enough lines to parse the positions
            # should be same format as input card
            n_atoms = len(prev_structure)
            positions_card = get_atomic_positions(
                pwo_lines[image_index:image_index + n_atoms + 1],
                n_atoms=n_atoms, cell=cell, alat=cell_alat)

            # convert to Atoms object
            symbols = [label_to_symbol(position[0]) for position in
                       positions_card]
            custom_labels = [position[0] for position in positions_card]
            positions = [position[1] for position in positions_card]
            constraint_flags = [position[2] for position in positions_card]
            structure = AtomsCustom(symbols=symbols, positions=positions, cell=cell,
                              pbc=True, custom_labels=custom_labels)
            structure.set_constraint(convert_constraint_flags(constraint_flags))

        # Extract calculation results
        # Energy
        energy = None
        for energy_index in indexes[_PW_TOTEN]:
            if image_index < energy_index < next_index:
                energy = float(
                    pwo_lines[energy_index].split()[-2]) * units['Ry']

        # Forces
        forces = None
        for force_index in indexes[_PW_FORCE]:
            if image_index < force_index < next_index:
                # Before QE 5.3 'negative rho' added 2 lines before forces
                # Use exact lines to stop before 'non-local' forces
                # in high verbosity
                if not pwo_lines[force_index + 2].strip():
                    force_index += 4
                else:
                    force_index += 2
                # assume contiguous
                forces = [
                    [float(x) for x in force_line.split()[-3:]] for force_line
                    in pwo_lines[force_index:force_index + len(structure)]]
                forces = np.array(forces) * units['Ry'] / units['Bohr']

        # Stress
        stress = None
        for stress_index in indexes[_PW_STRESS]:
            if image_index < stress_index < next_index:
                sxx, sxy, sxz = pwo_lines[stress_index + 1].split()[:3]
                _, syy, syz = pwo_lines[stress_index + 2].split()[:3]
                _, _, szz = pwo_lines[stress_index + 3].split()[:3]
                stress = np.array([sxx, syy, szz, syz, sxz, sxy], dtype=float)
                # sign convention is opposite of ase
                stress *= -1 * units['Ry'] / (units['Bohr'] ** 3)

        # Magmoms
        magmoms = None
        for magmoms_index in indexes[_PW_MAGMOM]:
            if image_index < magmoms_index < next_index:
                magmoms = [
                    float(mag_line.split()[-1]) for mag_line
                    in pwo_lines[magmoms_index + 1:
                                 magmoms_index + 1 + len(structure)]]

        # Dipole moment
        dipole = None
        if indexes[_PW_DIPOLE]:
            for dipole_index in indexes[_PW_DIPOLE]:
                if image_index < dipole_index < next_index:
                    _dipole = float(pwo_lines[dipole_index].split()[-2])

            for dipole_index in indexes[_PW_DIPOLE_DIRECTION]:
                if image_index < dipole_index < next_index:
                    _direction = pwo_lines[dipole_index].strip()
                    prefix = 'Computed dipole along edir('
                    _direction = _direction[len(prefix):]
                    _direction = int(_direction[0])

            if not _dipole:
                _dipole = 0

            dipole = np.eye(3)[_direction - 1] * _dipole * units['Debye']

        # Fermi level / highest occupied level
        efermi = None
        for fermi_index in indexes[_PW_FERMI]:
            if image_index < fermi_index < next_index:
                efermi = float(pwo_lines[fermi_index].split()[-2])

        if efermi is None:
            for ho_index in indexes[_PW_HIGHEST_OCCUPIED]:
                if image_index < ho_index < next_index:
                    efermi = float(pwo_lines[ho_index].split()[-1])

        if efermi is None:
            for holf_index in indexes[_PW_HIGHEST_OCCUPIED_LOWEST_FREE]:
                if image_index < holf_index < next_index:
                    efermi = float(pwo_lines[holf_index].split()[-2])

        # K-points
        ibzkpts = None
        weights = None
        kpoints_warning = "Number of k-points >= 100: " + \
                          "set verbosity='high' to print them."

        for kpts_index in indexes[_PW_KPTS]:
            nkpts = int(re.findall(r'\b\d+\b', pwo_lines[kpts_index])[0])
            kpts_index += 2

            if pwo_lines[kpts_index].strip() == kpoints_warning:
                continue

            # QE prints the k-points in units of 2*pi/alat
            # with alat defined as the length of the first
            # cell vector
            cell = structure.get_cell()
            alat = np.linalg.norm(cell[0])
            ibzkpts = []
            weights = []
            for i in range(nkpts):
                L = pwo_lines[kpts_index + i].split()
                weights.append(float(L[-1]))
                coord = np.array([L[-6], L[-5], L[-4].strip('),')],
                                 dtype=float)
                coord *= 2 * np.pi / alat
                coord = kpoint_convert(cell, ckpts_kv=coord)
                ibzkpts.append(coord)
            ibzkpts = np.array(ibzkpts)
            weights = np.array(weights)


        # Bands
        kpts = None
        kpoints_warning = "Number of k-points >= 100: " + \
                          "set verbosity='high' to print the bands."

        try:
            for bands_index in indexes[_PW_BANDS] + indexes[_PW_BANDSTRUCTURE]:
                if image_index < bands_index < next_index:
                    bands_index += 1
                    # skip over the lines with DFT+U occupation matrices
                    if 'enter write_ns' in pwo_lines[bands_index]:
                        while 'exit write_ns' not in pwo_lines[bands_index]:
                            bands_index += 1
                    bands_index += 1

                    if pwo_lines[bands_index].strip() == kpoints_warning:
                        continue

                    assert ibzkpts is not None
                    spin, bands, eigenvalues = 0, [], [[], []]

                    while True:
                        L = pwo_lines[bands_index].replace('-', ' -').split()
                        if len(L) == 0:
                            if len(bands) > 0:
                                eigenvalues[spin].append(bands)
                                bands = []
                        elif L == ['occupation', 'numbers']:
                            # Skip the lines with the occupation numbers
                            bands_index += len(eigenvalues[spin][0]) // 8 + 1
                        elif L[0] == 'k' and L[1].startswith('='):
                            pass
                        elif 'SPIN' in L:
                            if 'DOWN' in L:
                                spin += 1
                        else:
                            try:
                                bands.extend(map(float, L))
                            except ValueError:
                                break
                        bands_index += 1

                    if spin == 1:
                        assert len(eigenvalues[0]) == len(eigenvalues[1])
                    assert len(eigenvalues[0]) == len(ibzkpts), \
                        (np.shape(eigenvalues), len(ibzkpts))

                    kpts = []
                    for s in range(spin + 1):
                        for w, k, e in zip(weights, ibzkpts, eigenvalues[s]):
                            kpt = SinglePointKPoint(w, s, k, eps_n=e)
                            kpts.append(kpt)
        except AssertionError:
            print('Warning: bands were not read from pwo file.')


        # Put everything together
        #
        # In PW the forces are consistent with the "total energy"; that's why
        # its value must be assigned to free_energy.
        # PW doesn't compute the extrapolation of the energy to 0K smearing
        # the closer thing to this is again the total energy that contains
        # the correct (i.e. variational) form of the band energy is
        #   Eband = \int e N(e) de   for e<Ef , where N(e) is the DOS
        # This differs by the term (-TS)  from the sum of KS eigenvalues:
        #    Eks = \sum wg(n,k) et(n,k)
        # which is non variational. When a Fermi-Dirac function is used
        # for a given T, the variational energy is REALLY the free energy F,
        # and F = E - TS , with E = non variational energy.
        #
        calc = SinglePointDFTCalculator(structure, energy=energy,
                                        free_energy=energy,
                                        forces=forces, stress=stress,
                                        magmoms=magmoms, efermi=efermi,
                                        ibzkpts=ibzkpts, dipole=dipole)
        calc.kpts = kpts
        structure.calc = calc

        yield structure


# Runtime patching
ase.io.espresso.write_espresso_in = write_espresso_in_custom
ase.io.espresso.read_espresso_in = read_espresso_in_custom
ase.io.espresso.read_espresso_out = read_espresso_out_custom
