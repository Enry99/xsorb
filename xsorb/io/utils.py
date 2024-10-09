from ase.io import read

from xsorb import ase_custom

def ase_custom_read(filename, **kwargs):
    return ase_custom.Atoms_custom(read(filename=filename, **kwargs))
