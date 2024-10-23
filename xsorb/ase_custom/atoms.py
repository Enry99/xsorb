#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Enrico Pedretti
# Credit to original ASE code: https://wiki.fysik.dtu.dk/ase/

'''
Redefinition of the ase.Atoms class to:
- include custom labels for the atoms, obtainable with get_custom_labels()
- add constraints when summing two Atoms objects.
'''

import re

from ase import Atoms, Atom
from ase.constraints import FixAtoms, FixCartesian, FixScaled


class AtomsCustom(Atoms):
    """
    Custom class inheriting from ASE Atoms class,
    with the addition of custom_labels array,
    and the corresponding get_custom_labels() and set_custom_labels() methods.

    It also has the added functionality
    to include both constraints when summing two Atoms objects.
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
                 custom_labels=None):

        if hasattr(symbols, 'get_custom_labels'):
            #here symbols is a AtomsCustom:
            # copy constructor
            custom_labels = symbols.get_custom_labels()
        elif symbols is not None:
            #strip possible custom labels:
            for i, symbol in enumerate(symbols):
                if isinstance(symbol, str):
                    symbols[i] = re.sub(r'\d+', '', symbol) #remove numbers

        super().__init__(symbols=symbols,positions=positions,numbers=numbers,tags=tags,
                         momenta=momenta,masses=masses,magmoms=magmoms,charges=charges,
                         scaled_positions=scaled_positions,cell=cell, pbc=pbc, celldisp=celldisp,
                         constraint=constraint,calculator=calculator,info=info)

        if custom_labels is None:
            custom_labels = self.get_chemical_symbols()

        self.set_custom_labels(custom_labels)


    def set_custom_labels(self, custom_labels=None):
        """Set custom_labels."""

        if custom_labels is not None:
            labels = [f'{label:3}' for label in custom_labels]
        else: labels=None
        self.set_array('custom_labels', labels, str, ())


    def get_custom_labels(self):
        """Get array of custom_labels."""
        if 'custom_labels' in self.arrays:
            return [label.strip() for label in self.arrays['custom_labels']]
        else:
            return None


    def extend(self, other):
        """Extend atoms object by appending atoms from *other*.
        Also adds constraints from the second object.
        """

        if isinstance(other, Atom):
            other = self.__class__([other])
        othercp = other.copy()
        n1 = len(self)

        super().extend(other)

        #add constraints from the second
        for constr in othercp.constraints:
            if isinstance(constr, (FixCartesian,FixScaled,FixAtoms)):
                constr.a += n1
                self.constraints += [constr]
