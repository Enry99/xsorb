#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from ase.atoms import Atoms, default
from ase.io.espresso import *


import copy
import numpy as np
import ase.units as units
from ase.atom import Atom
from ase.cell import Cell
from typing import List
from ase.symbols import string2symbols

#this class stores the custom labels of the atoms as a list, obtainable with get_custom_labels()
class Atoms_custom(Atoms):
    """Atoms object.

    The Atoms object can represent an isolated molecule, or a
    periodically repeated structure.  It has a unit cell and
    there may be periodic boundary conditions along any of the three
    unit cell axes.
    Information about the atoms (atomic numbers and position) is
    stored in ndarrays.  Optionally, there can be information about
    tags, momenta, masses, magnetic moments and charges.

    In order to calculate energies, forces and stresses, a calculator
    object has to attached to the atoms object.

    Parameters:

    symbols: str (formula) or list of str
        Can be a string formula, a list of symbols or a list of
        Atom objects.  Examples: 'H2O', 'COPt12', ['H', 'H', 'O'],
        [Atom('Ne', (x, y, z)), ...].
    positions: list of xyz-positions
        Atomic positions.  Anything that can be converted to an
        ndarray of shape (n, 3) will do: [(x1,y1,z1), (x2,y2,z2),
        ...].
    scaled_positions: list of scaled-positions
        Like positions, but given in units of the unit cell.
        Can not be set at the same time as positions.
    numbers: list of int
        Atomic numbers (use only one of symbols/numbers).
    tags: list of int
        Special purpose tags.
    momenta: list of xyz-momenta
        Momenta for all atoms.
    masses: list of float
        Atomic masses in atomic units.
    magmoms: list of float or list of xyz-values
        Magnetic moments.  Can be either a single value for each atom
        for collinear calculations or three numbers for each atom for
        non-collinear calculations.
    charges: list of float
        Initial atomic charges.
    cell: 3x3 matrix or length 3 or 6 vector
        Unit cell vectors.  Can also be given as just three
        numbers for orthorhombic cells, or 6 numbers, where
        first three are lengths of unit cell vectors, and the
        other three are angles between them (in degrees), in following order:
        [len(a), len(b), len(c), angle(b,c), angle(a,c), angle(a,b)].
        First vector will lie in x-direction, second in xy-plane,
        and the third one in z-positive subspace.
        Default value: [0, 0, 0].
    celldisp: Vector
        Unit cell displacement vector. To visualize a displaced cell
        around the center of mass of a Systems of atoms. Default value
        = (0,0,0)
    pbc: one or three bool
        Periodic boundary conditions flags.  Examples: True,
        False, 0, 1, (1, 1, 0), (True, False, False).  Default
        value: False.
    constraint: constraint object(s)
        Used for applying one or more constraints during structure
        optimization.
    calculator: calculator object
        Used to attach a calculator for calculating energies and atomic
        forces.
    info: dict of key-value pairs
        Dictionary of key-value pairs with additional information
        about the system.  The following keys may be used by ase:

          - spacegroup: Spacegroup instance
          - unit_cell: 'conventional' | 'primitive' | int | 3 ints
          - adsorbate_info: Information about special adsorption sites

        Items in the info attribute survives copy and slicing and can
        be stored in and retrieved from trajectory files given that the
        key is a string, the value is JSON-compatible and, if the value is a
        user-defined object, its base class is importable.  One should
        not make any assumptions about the existence of keys.

    Examples:

    These three are equivalent:

    >>> d = 1.104  # N2 bondlength
    >>> a = Atoms('N2', [(0, 0, 0), (0, 0, d)])
    >>> a = Atoms(numbers=[7, 7], positions=[(0, 0, 0), (0, 0, d)])
    >>> a = Atoms([Atom('N', (0, 0, 0)), Atom('N', (0, 0, d))])

    FCC gold:

    >>> a = 4.05  # Gold lattice constant
    >>> b = a / 2
    >>> fcc = Atoms('Au',
    ...             cell=[(0, b, b), (b, 0, b), (b, b, 0)],
    ...             pbc=True)

    Hydrogen wire:

    >>> d = 0.9  # H-H distance
    >>> h = Atoms('H', positions=[(0, 0, 0)],
    ...           cell=(d, 0, 0),
    ...           pbc=(1, 0, 0))
    """

    def __init__(self, symbols=None,
                 positions=None, numbers=None,
                 tags=None, momenta=None, masses=None,
                 magmoms=None, charges=None,
                 scaled_positions=None,
                 cell=None, pbc=None, celldisp=None,
                 constraint=None,
                 calculator=None,
                 info=None,
                 velocities=None,
                 custom_labels=None):

        self._cellobj = Cell.new()
        self._pbc = np.zeros(3, bool)

        atoms = None

        if hasattr(symbols, 'get_positions'):
            atoms = symbols
            symbols = None
        elif (isinstance(symbols, (list, tuple)) and
              len(symbols) > 0 and isinstance(symbols[0], Atom)):
            # Get data from a list or tuple of Atom objects:
            data = [[atom.get_raw(name) for atom in symbols]
                    for name in
                    ['position', 'number', 'tag', 'momentum',
                     'mass', 'magmom', 'charge']]
            atoms = self.__class__(None, *data)
            symbols = None

        if atoms is not None:
            # Get data from another Atoms object:
            if scaled_positions is not None:
                raise NotImplementedError
            if symbols is None and numbers is None:
                numbers = atoms.get_atomic_numbers()
            if positions is None:
                positions = atoms.get_positions()
            if tags is None and atoms.has('tags'):
                tags = atoms.get_tags()
            if momenta is None and atoms.has('momenta'):
                momenta = atoms.get_momenta()
            if magmoms is None and atoms.has('initial_magmoms'):
                magmoms = atoms.get_initial_magnetic_moments()
            if masses is None and atoms.has('masses'):
                masses = atoms.get_masses()
            if charges is None and atoms.has('initial_charges'):
                charges = atoms.get_initial_charges()
            if cell is None:
                cell = atoms.get_cell()
            if celldisp is None:
                celldisp = atoms.get_celldisp()
            if pbc is None:
                pbc = atoms.get_pbc()
            if constraint is None:
                constraint = [c.copy() for c in atoms.constraints]
            if calculator is None:
                calculator = atoms.calc
            if info is None:
                info = copy.deepcopy(atoms.info)
            if custom_labels is None and atoms.has('custom_labels'):
                custom_labels = atoms.get_custom_labels()

        self.arrays = {}

        if symbols is None:
            if numbers is None:
                if positions is not None:
                    natoms = len(positions)
                elif scaled_positions is not None:
                    natoms = len(scaled_positions)
                else:
                    natoms = 0
                numbers = np.zeros(natoms, int)
            self.new_array('numbers', numbers, int)
        else:
            if numbers is not None:
                raise TypeError(
                    'Use only one of "symbols" and "numbers".')
            else:
                self.new_array('numbers', self.symbols2numbers(symbols), int)

        if self.numbers.ndim != 1:
            raise ValueError('"numbers" must be 1-dimensional.')

        if cell is None:
            cell = np.zeros((3, 3))
        self.set_cell(cell)

        if celldisp is None:
            celldisp = np.zeros(shape=(3, 1))
        self.set_celldisp(celldisp)

        if positions is None:
            if scaled_positions is None:
                positions = np.zeros((len(self.arrays['numbers']), 3))
            else:
                assert self.cell.rank == 3
                positions = np.dot(scaled_positions, self.cell)
        else:
            if scaled_positions is not None:
                raise TypeError(
                    'Use only one of "symbols" and "numbers".')
        self.new_array('positions', positions, float, (3,))

        self.set_constraint(constraint)
        self.set_tags(default(tags, 0))
        self.set_masses(default(masses, None))
        self.set_initial_magnetic_moments(default(magmoms, 0.0))
        self.set_initial_charges(default(charges, 0.0))
        if pbc is None:
            pbc = False
        self.set_pbc(pbc)
        self.set_momenta(default(momenta, (0.0, 0.0, 0.0)),
                         apply_constraint=False)

        if velocities is not None:
            if momenta is None:
                self.set_velocities(velocities)
            else:
                raise TypeError(
                    'Use only one of "momenta" and "velocities".')

        if info is None:
            self.info = {}
        else:
            self.info = dict(info)

        self.calc = calculator

        self.set_custom_labels(custom_labels)

    def set_custom_labels(self, custom_labels=None):
        """Set custom_labels."""
        if custom_labels is None:
            self.set_array('custom_labels', None)
        self.set_array('custom_labels', custom_labels, str, ())

    def get_custom_labels(self):
        """Get array of custom_labels."""
        if 'custom_labels' in self.arrays:
            return self.arrays['custom_labels'].copy()
        else:
            return None

    def extend(self, other):
        """Extend atoms object by appending atoms from *other*."""
        if isinstance(other, Atom):
            other = self.__class__([other])

        n1 = len(self)
        n2 = len(other)

        for name, a1 in self.arrays.items():
            a = np.zeros((n1 + n2,) + a1.shape[1:], a1.dtype)
            a[:n1] = a1
            if name == 'masses':
                a2 = other.get_masses()
            else:
                a2 = other.arrays.get(name)
            if a2 is not None:
                a[n1:] = a2
            self.arrays[name] = a

        for name, a2 in other.arrays.items():
            if name in self.arrays:
                continue
            a = np.empty((n1 + n2,) + a2.shape[1:], a2.dtype)
            a[n1:] = a2
            if name == 'masses':
                a[:n1] = self.get_masses()[:n1]
            else:
                a[:n1] = 0

            self.set_array(name, a)

        #add constraints from the second
        othercp = other.copy()
        for constr in othercp.constraints:
            if isinstance(constr, FixCartesian):  
                constr.a += n1         
                self.constraints += [constr]

    def symbols2numbers(self, symbols) -> List[int]:
        if isinstance(symbols, str):
            symbols = string2symbols(symbols)
        numbers = []
        for s in symbols:
            if isinstance(s, str):
                if(len(s)>2): s = s[:2]
                if(len(s)>1):
                    if(s[1].isdigit()):
                        s = s[0]            
                numbers.append(atomic_numbers[s])
            else:
                numbers.append(int(s))
        return numbers
#################################################

#simply reads labels for the atoms and stores them in an interal list.
@iofunction('rU')
def read_espresso_in_custom(fileobj):
    """Parse a Quantum ESPRESSO input files, '.in', '.pwi'.

    ESPRESSO inputs are generally a fortran-namelist format with custom
    blocks of data. The namelist is parsed as a dict and an atoms object
    is constructed from the included information.

    Parameters
    ----------
    fileobj : file | str
        A file-like object that supports line iteration with the contents
        of the input file, or a filename.

    Returns
    -------
    atoms : Atoms
        Structure defined in the input file.

    Raises
    ------
    KeyError
        Raised for missing keys that are required to process the file
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
        cell, cell_alat = get_cell_parameters(card_lines, alat=alat)
    else:
        alat, cell = ibrav_to_cell(data['system'])

    # species_info holds some info for each element
    species_card = get_atomic_species(card_lines, n_species=data['system']['ntyp'])
    species_info = {}
    for ispec, (label, weight, pseudo) in enumerate(species_card):
        symbol = label_to_symbol(label)
        valence = get_valence_electrons(symbol, data, pseudo)

        # starting_magnetization is in fractions of valence electrons
        magnet_key = "starting_magnetization({0})".format(ispec + 1)
        magmom = valence * data["system"].get(magnet_key, 0.0)
        species_info[label] = {"weight": weight, "pseudo": pseudo,
                                "valence": valence, "magmom": magmom}
    

    positions_card = get_atomic_positions(
        card_lines, n_atoms=data['system']['nat'], cell=cell, alat=alat)

    symbols = [label_to_symbol(position[0]) for position in positions_card]
    custom_labels = [position[0] for position in positions_card]
    positions = [position[1] for position in positions_card]
    magmoms = [species_info[label]["magmom"] for label in custom_labels]
    

    # TODO: put more info into the atoms object
    # e.g magmom, forces.
    atoms = Atoms_custom(symbols=symbols, positions=positions, cell=cell, pbc=True,
                  magmoms=magmoms, custom_labels=custom_labels)

    return atoms

#################################################

#This simply passes through the pseudopotentials, leaving the order unchanged, and keeps the labels of the atoms
def write_espresso_in_custom(fd, atoms, input_data=None, pseudopotentials=None,
                      kspacing=None, kpts=None, koffset=(0, 0, 0),
                      crystal_coordinates=False, **kwargs):
    """
    Create an input file for pw.x.

    Use set_initial_magnetic_moments to turn on spin, if ispin is set to 2
    with no magnetic moments, they will all be set to 0.0. Magnetic moments
    will be converted to the QE units (fraction of valence electrons) using
    any pseudopotential files found, or a best guess for the number of
    valence electrons.

    Units are not converted for any other input data, so use Quantum ESPRESSO
    units (Usually Ry or atomic units).

    Keys with a dimension (e.g. Hubbard_U(1)) will be incorporated as-is
    so the `i` should be made to match the output.

    Implemented features:

    - Conversion of :class:`ase.constraints.FixAtoms` and
                    :class:`ase.constraints.FixCartesian`.
    - `starting_magnetization` derived from the `mgmoms` and pseudopotentials
      (searches default paths for pseudo files.)
    - Automatic assignment of options to their correct sections.
    - Interpretation of ibrav (cell must exactly match the vectors defined
      in the QE docs).

    Not implemented:

    - Lists of k-points
    - Other constraints
    - Hubbard parameters
    - Validation of the argument types for input
    - Validation of required options
    - Reorientation for ibrav settings
    - Noncollinear magnetism

    Parameters
    ----------
    fd: file
        A file like object to write the input file to.
    atoms: Atoms
        A single atomistic configuration to write to `fd`.
    input_data: dict
        A flat or nested dictionary with input parameters for pw.x
    pseudopotentials: dict
        A filename for each atomic species, e.g.
        {'O': 'O.pbe-rrkjus.UPF', 'H': 'H.pbe-rrkjus.UPF'}.
        A dummy name will be used if none are given.
    kspacing: float
        Generate a grid of k-points with this as the minimum distance,
        in A^-1 between them in reciprocal space. If set to None, kpts
        will be used instead.
    kpts: (int, int, int) or dict
        If kpts is a tuple (or list) of 3 integers, it is interpreted
        as the dimensions of a Monkhorst-Pack grid.
        If ``kpts`` is set to ``None``, only the Γ-point will be included
        and QE will use routines optimized for Γ-point-only calculations.
        Compared to Γ-point-only calculations without this optimization
        (i.e. with ``kpts=(1, 1, 1)``), the memory and CPU requirements
        are typically reduced by half.
        If kpts is a dict, it will either be interpreted as a path
        in the Brillouin zone (*) if it contains the 'path' keyword,
        otherwise it is converted to a Monkhorst-Pack grid (**).
        (*) see ase.dft.kpoints.bandpath
        (**) see ase.calculators.calculator.kpts2sizeandoffsets
    koffset: (int, int, int)
        Offset of kpoints in each direction. Must be 0 (no offset) or
        1 (half grid offset). Setting to True is equivalent to (1, 1, 1).
    crystal_coordinates: bool
        Whether the atomic positions should be written to the QE input file in
        absolute (False, default) or relative (crystal) coordinates (True).

    """


    # Convert to a namelist to make working with parameters much easier
    # Note that the name ``input_data`` is chosen to prevent clash with
    # ``parameters`` in Calculator objects
    input_parameters = construct_namelist_custom(input_data, **kwargs)

    # Convert ase constraints to QE constraints
    # Nx3 array of force multipliers matches what QE uses
    # Do this early so it is available when constructing the atoms card
    constraint_mask = np.ones((len(atoms), 3), dtype='int')
    for constraint in atoms.constraints:
        if isinstance(constraint, FixAtoms):
            constraint_mask[constraint.index] = 0
        elif isinstance(constraint, FixCartesian):
            constraint_mask[constraint.a] = constraint.mask
        else:
            warnings.warn('Ignored unknown constraint {}'.format(constraint))

    # Species info holds the information on the pseudopotential and
    # associated for each element
    #species_info = {}
    #for label, pseudo in pseudopotentials.items():
    #    species_info[label] = {'pseudo': pseudo, 'valence': None}

    # Convert atoms into species.
    # Each different magnetic moment needs to be a separate type even with
    # the same pseudopotential (e.g. an up and a down for AFM).
    # if any magmom are > 0 or nspin == 2 then use species labels.
    # Rememeber: magnetisation uses 1 based indexes
    atomic_species = OrderedDict()
    atomic_species_str = []
    atomic_positions_str = []

    
    if False: #do not force spin on, simply passtrhough
        nspin = input_parameters['system'].get('nspin', 1)  # 1 is the default
        if any(atoms.get_initial_magnetic_moments()):
            if nspin == 1:
                # Force spin on
                input_parameters['system']['nspin'] = 2
                nspin = 2

    if False: #nspin == 2:
        # Spin on
        for atom, magmom in zip(atoms, atoms.get_initial_magnetic_moments()):
            if (atom.symbol, magmom) not in atomic_species:
                # spin as fraction of valence
                fspin = float(magmom) / species_info[atom.symbol]['valence']
                # Index in the atomic species list
                sidx = len(atomic_species) + 1
                # Index for that atom type; no index for first one
                tidx = sum(atom.symbol == x[0] for x in atomic_species) or ' '
                atomic_species[(atom.symbol, magmom)] = (sidx, tidx)
                # Add magnetization to the input file
                mag_str = 'starting_magnetization({0})'.format(sidx)
                input_parameters['system'][mag_str] = fspin
                atomic_species_str.append(
                    '{species}{tidx} {mass} {pseudo}\n'.format(
                        species=atom.symbol, tidx=tidx, mass=atom.mass,
                        pseudo=species_info[atom.symbol]['pseudo']))
            # lookup tidx to append to name
            sidx, tidx = atomic_species[(atom.symbol, magmom)]

            # only inclued mask if something is fixed
            if not all(constraint_mask[atom.index]):
                mask = ' {mask[0]} {mask[1]} {mask[2]}'.format(
                    mask=constraint_mask[atom.index])
            else:
                mask = ''

            # construct line for atomic positions
            atomic_positions_str.append(
                '{atom.symbol}{tidx} '
                '{atom.x:.10f} {atom.y:.10f} {atom.z:.10f}'
                '{mask}\n'.format(atom=atom, tidx=tidx, mask=mask))

    else:
        # Do nothing about magnetisation

        for label, pseudo in pseudopotentials.items():
            atomic_species_str.append(
                    '{label} {mass} {pseudo}\n'.format(
                        label=label, 
                        mass=Atom(label_to_symbol(label)).mass,
                        pseudo=pseudo))

        for atom, label in zip(atoms, atoms.get_custom_labels()):

            # only inclued mask if something is fixed
            if not all(constraint_mask[atom.index]):
                mask = ' {mask[0]} {mask[1]} {mask[2]}'.format(
                    mask=constraint_mask[atom.index])
            else:
                mask = ''

            if crystal_coordinates:
                coords = [atom.a, atom.b, atom.c]
            else:
                coords = atom.position
            atomic_positions_str.append(
                '{label} '
                '{coords[0]:.10f} {coords[1]:.10f} {coords[2]:.10f} '
                '{mask}\n'.format(label=label, coords=coords, mask=mask))

    # Add computed parameters
    # different magnetisms means different types
    input_parameters['system']['ntyp'] = len(pseudopotentials)
    input_parameters['system']['nat'] = len(atoms)

    # Use cell as given or fit to a specific ibrav
    if 'ibrav' in input_parameters['system']:
        ibrav = input_parameters['system']['ibrav']
        if ibrav != 0:
            celldm = cell_to_ibrav(atoms.cell, ibrav)
            regen_cell = ibrav_to_cell(celldm)[1]
            if not np.allclose(atoms.cell, regen_cell):
                warnings.warn('Input cell does not match requested ibrav'
                              '{} != {}'.format(regen_cell, atoms.cell))
            input_parameters['system'].update(celldm)
    else:
        # Just use standard cell block
        input_parameters['system']['ibrav'] = 0

    # Construct input file into this
    pwi = []

    # Assume sections are ordered (taken care of in namelist construction)
    # and that repr converts to a QE readable representation (except bools)
    for section in input_parameters:
        pwi.append('&{0}\n'.format(section.upper()))
        for key, value in input_parameters[section].items():
            if value is True:
                pwi.append('   {0:16} = .true.\n'.format(key))
            elif value is False:
                pwi.append('   {0:16} = .false.\n'.format(key))
            else:
                # repr format to get quotes around strings
                pwi.append('   {0:16} = {1!r:}\n'.format(key, value))
        pwi.append('/\n')  # terminate section
    pwi.append('\n')

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
        pwi.append('%s\n' % len(kgrid))
        for k in kgrid:
            pwi.append('{k[0]:.14f} {k[1]:.14f} {k[2]:.14f} 0\n'.format(k=k))
        pwi.append('\n')
    elif isinstance(kgrid, str) and (kgrid == "gamma"):
        pwi.append('K_POINTS gamma\n')
        pwi.append('\n')
    else:
        pwi.append('K_POINTS automatic\n')
        pwi.append('{0[0]} {0[1]} {0[2]}  {1[0]:d} {1[1]:d} {1[2]:d}\n'
                   ''.format(kgrid, koffset))
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

#bugfix for OrderedDict mutated
def construct_namelist_custom(parameters=None, warn=False, **kwargs):
    """
    Construct an ordered Namelist containing all the parameters given (as
    a dictionary or kwargs). Keys will be inserted into their appropriate
    section in the namelist and the dictionary may contain flat and nested
    structures. Any kwargs that match input keys will be incorporated into
    their correct section. All matches are case-insensitive, and returned
    Namelist object is a case-insensitive dict.

    If a key is not known to ase, but in a section within `parameters`,
    it will be assumed that it was put there on purpose and included
    in the output namelist. Anything not in a section will be ignored (set
    `warn` to True to see ignored keys).

    Keys with a dimension (e.g. Hubbard_U(1)) will be incorporated as-is
    so the `i` should be made to match the output.

    The priority of the keys is:
        kwargs[key] > parameters[key] > parameters[section][key]
    Only the highest priority item will be included.

    Parameters
    ----------
    parameters: dict
        Flat or nested set of input parameters.
    warn: bool
        Enable warnings for unused keys.

    Returns
    -------
    input_namelist: Namelist
        pw.x compatible namelist of input parameters.

    """
    # Convert everything to Namelist early to make case-insensitive
    if parameters is None:
        parameters = Namelist()
    else:
        # Maximum one level of nested dict
        # Don't modify in place
        parameters_namelist = Namelist()
        for key, value in parameters.items():
            if isinstance(value, dict):
                parameters_namelist[key] = Namelist(value)
            else:
                parameters_namelist[key] = value
        parameters = parameters_namelist

    # Just a dict
    kwargs = Namelist(kwargs)

    # Final parameter set
    input_namelist = Namelist()

    # Collect
    for section in KEYS:
        sec_list = Namelist()
        for key in KEYS[section]:
            # Check all three separately and pop them all so that
            # we can check for missing values later
            # This first block checks only keys not in the format (i)
            if key in parameters.get(section, {}):  #if given in a dictionary inside sections
                sec_list[key] = parameters[section].pop(key)
            if key in parameters:  #if passed directly as a named parameter to the function, e.g. pseudopotentials=...
                sec_list[key] = parameters.pop(key)
            if key in kwargs: #if passed as a kwarg
                sec_list[key] = kwargs.pop(key)

            # Check if there is a key(i) version (no extra parsing)
            cp_parameters_section = parameters.get(section, {}).copy()  #NOTE: bugfix to ase source, since we cannot pop on the Namelist that we are looping on
            #NOTE: It will probably be fixed in next ase release (3.23)
            for arg_key in cp_parameters_section:
                if arg_key.split('(')[0].strip().lower() == key.lower():
                    sec_list[arg_key] = parameters[section].pop(arg_key)                 
            cp_parameters = parameters.copy()
            for arg_key in cp_parameters:
                if arg_key.split('(')[0].strip().lower() == key.lower():
                    sec_list[arg_key] = parameters.pop(arg_key)
            cp_kwargs = kwargs.copy()
            for arg_key in cp_kwargs:
                if arg_key.split('(')[0].strip().lower() == key.lower():
                    sec_list[arg_key] = kwargs.pop(arg_key)

        # Add to output
        input_namelist[section] = sec_list

    unused_keys = list(kwargs)
    # pass anything else already in a section
    for key, value in parameters.items():
        if key in KEYS and isinstance(value, dict):
            input_namelist[key].update(value)
        elif isinstance(value, dict):
            unused_keys.extend(list(value))
        else:
            unused_keys.append(key)

    if warn and unused_keys:
        warnings.warn('Unused keys: {}'.format(', '.join(unused_keys)))

    return input_namelist



#################################################

units = create_units('2006')
# Section identifiers
_PW_START = 'Program PWSCF'
_PW_END = 'End of self-consistent calculation'
_PW_CELL = 'CELL_PARAMETERS'
_PW_POS = 'ATOMIC_POSITIONS'
_PW_MAGMOM = 'Magnetic moment per site'
_PW_FORCE = 'Forces acting on atoms'
_PW_TOTEN = '!    total energy'
_PW_STRESS = 'total   stress'
_PW_FERMI = 'the Fermi energy is'
_PW_HIGHEST_OCCUPIED = 'highest occupied level'
_PW_HIGHEST_OCCUPIED_LOWEST_FREE = 'highest occupied, lowest unoccupied level'
_PW_KPTS = 'number of k points='
_PW_BANDS = _PW_END
_PW_BANDSTRUCTURE = 'End of band structure calculation'

def read_espresso_out_custom(fileobj, index=-1, results_required=True):
    """Reads Quantum ESPRESSO output files.

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
            if any([config_index < results_index < config_index_next
                    for results_index in results_indexes]):
                results_config_indexes.append(config_index)

        # slice from the subset
        image_indexes = results_config_indexes[index]
    else:
        image_indexes = all_config_indexes[index]

    # Extract initialisation information each time PWSCF starts
    # to add to subsequent configurations. Use None so slices know
    # when to fill in the blanks.
    pwscf_start_info = dict((idx, None) for idx in indexes[_PW_START])

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
            structure = Atoms_custom(symbols=symbols, positions=positions, cell=cell,
                              pbc=True, custom_labels=custom_labels)

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
            nkpts = int(pwo_lines[kpts_index].split()[4])
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

        for bands_index in indexes[_PW_BANDS] + indexes[_PW_BANDSTRUCTURE]:
            if image_index < bands_index < next_index:
                bands_index += 2

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

        # Put everything together
        #
        # I have added free_energy.  Can and should we distinguish
        # energy and free_energy?  --askhl
        calc = SinglePointDFTCalculator(structure, energy=energy,
                                        free_energy=energy,
                                        forces=forces, stress=stress,
                                        magmoms=magmoms, efermi=efermi,
                                        ibzkpts=ibzkpts)
        calc.kpts = kpts
        structure.calc = calc

        yield structure


def parse_pwo_start_custom(lines, index=0):
    """Parse Quantum ESPRESSO calculation info from lines,
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
    # TODO: extend with extra DFT info?

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
    info['atoms'] = Atoms_custom(symbols=info['symbols'],
                          positions=info['positions'],
                          cell=info['cell'], pbc=True, custom_labels=custom_labels)

    return info


#################################################
from ase.io.extxyz import XYZError, Calculator, SinglePointCalculator, \
        key_val_str_to_dict, parse_properties, paropen, output_column_format,  \
        all_properties, per_atom_properties, per_config_properties

def _read_xyz_frame_custom(lines, natoms, properties_parser=key_val_str_to_dict,
                    nvec=0):
    # comment line
    line = next(lines).strip()
    if nvec > 0:
        info = {'comment': line}
    else:
        info = properties_parser(line) if line else {}

    pbc = None
    if 'pbc' in info:
        pbc = info['pbc']
        del info['pbc']
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
    for ln in range(natoms):
        try:
            line = next(lines)
        except StopIteration:
            raise XYZError('ase.io.extxyz: Frame has {} atoms, expected {}'
                           .format(len(data), natoms))
        vals = line.split()
        row = tuple([conv(val) for conv, val in zip(convs, vals)])
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
                raise XYZError('Expected cell vector, got {}'.format(entry[0]))

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

    symbols = None
    if 'symbols' in arrays:
        symbols = [s.capitalize() for s in arrays['symbols']]
        del arrays['symbols']

    numbers = None
    duplicate_numbers = None
    if 'numbers' in arrays:
        if symbols is None:
            numbers = arrays['numbers']
        else:
            duplicate_numbers = arrays['numbers']
        del arrays['numbers']

    charges = None
    if 'charges' in arrays:
        charges = arrays['charges']
        del arrays['charges']

    positions = None
    if 'positions' in arrays:
        positions = arrays['positions']
        del arrays['positions']

    atoms = Atoms_custom(symbols=symbols,
                  positions=positions,
                  numbers=numbers,
                  charges=charges,
                  cell=cell,
                  pbc=pbc,
                  info=info,
                  custom_labels=symbols)

    # Read and set constraints
    if 'move_mask' in arrays:
        if properties['move_mask'][1] == 3:
            cons = []
            for a in range(natoms):
                cons.append(FixCartesian(a, mask=~arrays['move_mask'][a, :]))
            atoms.set_constraint(cons)
        elif properties['move_mask'][1] == 1:
            atoms.set_constraint(FixAtoms(mask=~arrays['move_mask']))
        else:
            raise XYZError('Not implemented constraint')
        del arrays['move_mask']

    for name, array in arrays.items():
        atoms.new_array(name, array)

    if duplicate_numbers is not None:
        atoms.set_atomic_numbers(duplicate_numbers)

    # Load results of previous calculations into SinglePointCalculator
    results = {}
    for key in list(atoms.info.keys()):
        if key in per_config_properties:
            results[key] = atoms.info[key]
            # special case for stress- convert to Voigt 6-element form
            if key == 'stress' and results[key].shape == (3, 3):
                stress = results[key]
                stress = np.array([stress[0, 0],
                                   stress[1, 1],
                                   stress[2, 2],
                                   stress[1, 2],
                                   stress[0, 2],
                                   stress[0, 1]])
                results[key] = stress
    for key in list(atoms.arrays.keys()):
        if (key in per_atom_properties and len(value.shape) >= 1
            and value.shape[0] == len(atoms)):
            results[key] = atoms.arrays[key]
    if results != {}:
        calculator = SinglePointCalculator(atoms, **results)
        atoms.calc = calculator
    return atoms

def write_xyz_custom(fileobj, images, comment='', columns=None,
              write_info=True,
              write_results=True, plain=False, vec_cell=False,
              append=False, custom_labels_as_symbols = True):
    """
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
    if isinstance(fileobj, str):
        mode = 'w'
        if append:
            mode = 'a'
        fileobj = paropen(fileobj, mode)

    if hasattr(images, 'get_positions'):
        images = [images]

    for atoms in images:
        natoms = len(atoms)

        if columns is None:
            fr_cols = None
        else:
            fr_cols = columns[:]

        if fr_cols is None:
            fr_cols = (['symbols', 'positions']
                       + [key for key in atoms.arrays.keys() if
                          key not in ['symbols', 'positions', 'numbers',
                                      'species', 'pos']])

        if vec_cell:
            plain = True

        if plain:
            fr_cols = ['symbols', 'positions']
            write_info = False
            write_results = False

        per_frame_results = {}
        per_atom_results = {}
        if write_results:
            calculator = atoms.calc
            if (calculator is not None
                    and isinstance(calculator, Calculator)):
                for key in all_properties:
                    value = calculator.results.get(key, None)
                    if value is None:
                        # skip missing calculator results
                        continue
                    if (key in per_atom_properties and len(value.shape) >= 1
                        and value.shape[0] == len(atoms)):
                        # per-atom quantities (forces, energies, stresses)
                        per_atom_results[key] = value
                    elif key in per_config_properties:
                        # per-frame quantities (energy, stress)
                        # special case for stress, which should be converted
                        # to 3x3 matrices before writing
                        if key == 'stress':
                            xx, yy, zz, yz, xz, xy = value
                            value = np.array(
                                [(xx, xy, xz), (xy, yy, yz), (xz, yz, zz)])
                        per_frame_results[key] = value

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
            symbols = atoms.get_chemical_symbols()

        if natoms > 0 and not isinstance(symbols[0], str):
            raise ValueError('First column must be symbols-like')

        # Check second column "looks like" atomic positions
        pos = atoms.arrays[fr_cols[1]]
        if pos.shape != (natoms, 3) or pos.dtype.kind != 'f':
            raise ValueError('Second column must be position-like')

        # if vec_cell add cell information as pseudo-atoms
        if vec_cell:
            pbc = list(atoms.get_pbc())
            cell = atoms.get_cell()

            if True in pbc:
                nPBC = 0
                for i, b in enumerate(pbc):
                    if b:
                        nPBC += 1
                        symbols.append('VEC' + str(nPBC))
                        pos = np.vstack((pos, cell[i]))
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
                        cnstr[idx] = False
                elif isinstance(c0, FixCartesian):
                    masks = np.ones((natoms, 3), dtype=bool)
                    for i in range(len(cnstr)):
                        idx = cnstr[i].a
                        masks[idx] = cnstr[i].mask
                    cnstr = masks
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
                raise ValueError('Missing array "%s"' % column)

        #replace symbols with custom_labels
        if 'custom_labels' in arrays and custom_labels_as_symbols:
            arrays['symbols'] = arrays['custom_labels']
            del arrays['custom_labels']
            fr_cols.remove('custom_labels')

        if write_results:
            for key in per_atom_results:
                if key not in fr_cols:
                    fr_cols += [key]
                else:
                    warnings.warn('write_xyz() overwriting array "{0}" present '
                                  'in atoms.arrays with stored results '
                                  'from calculator'.format(key))
            arrays.update(per_atom_results)

        comm, ncols, dtype, fmt = output_column_format(atoms,
                                                       fr_cols,
                                                       arrays,
                                                       write_info,
                                                       per_frame_results)

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
        fileobj.write('%s\n' % comm)
        for i in range(natoms):
            fileobj.write(fmt % tuple(data[i]))



# Runtime patching
import ase.io.espresso
import ase.io.extxyz
ase.io.espresso.write_espresso_in = write_espresso_in_custom
ase.io.espresso.read_espresso_in = read_espresso_in_custom
ase.io.espresso.read_espresso_out = read_espresso_out_custom
ase.io.extxyz._read_xyz_frame = _read_xyz_frame_custom
ase.io.extxyz.write_xyz = write_xyz_custom
