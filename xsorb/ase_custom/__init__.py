'''
Mods to ase library to handle custom labels, and for rendering with POVRAY
'''

# Runtime patch for read/write
from xsorb.ase_custom import atoms, espresso, vaps, xyz

# Runtime patch for rendering is not called here but only explicitly
# with rom xsorb.ase_custom import povray when necessary
