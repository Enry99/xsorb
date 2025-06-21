#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Enrico Pedretti
# Credit to original ASE code: https://wiki.fysik.dtu.dk/ase/

'''
Redefinition of the ase.Atoms class to:
- include custom_labels property to atoms, handled by the tags property.
- add constraints when summing two Atoms objects.
'''

from ase import Atoms, Atom
from ase.constraints import FixAtoms, FixCartesian, FixScaled


class AtomsCustom(Atoms):
    """
    Custom class inheriting from ASE Atoms class,
    with custom labels included as tags.
    DO NOT USE TAGS SINCE THEY SHOULD BE RESERVED FOR CUSTOM LABELS.

    It also has the added functionality
    to include both constraints when summing two Atoms objects.
    """

    @property
    def custom_labels(self):
        """Get custom_labels."""

        symbols = self.get_chemical_symbols()
        tags = self.get_tags()

        if tags.any(): # if not all zero (default)
            custom_labels = []
            for symbol, tag in zip(self.get_chemical_symbols(), self.get_tags()):
                if tag != -1:
                    custom_labels.append(f"{symbol}{tag}")
                else:
                    custom_labels.append(symbol)
        else:
            custom_labels = symbols

        return custom_labels


    @custom_labels.setter
    def custom_labels(self, ids=None):
        """Set custom_labels."""
        if isinstance(ids, list) and all(isinstance(i, str) for i in ids):
            # if list of strings, extract numbers
            ids = [extract_number_from_string(s, self.symbols[i]) for i, s in enumerate(ids)]
            print(ids)
        self.set_tags(ids)


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
                constr.index[0] += n1
                self.constraints += [constr]


def extract_number_from_string(string, symbol):
    """
    Extract the number from a string, given the symbol.
    The number is the last part of the string, after the symbol.

    Parameters
    ----------
    string : str
        The string to extract the number from.
    symbol : str
        The symbol of the element.

    Returns
    -------
    number : int
        The number extracted from the string, or -1 if no number is found.
    """
    # extract the number from the end of the string
    number = string.split(symbol)[-1]
    if number:
        number = int(number)
    else:
        number = -1

    return number
