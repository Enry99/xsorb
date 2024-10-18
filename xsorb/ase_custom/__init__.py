'''
Mods to ase library to handle custom labels, and for rendering with POVRAY.
Need to be monkey patched at runtime by importin this module.
'''

# Runtime patch for read/write
import xsorb.ase_custom.espresso
import xsorb.ase_custom.vasp
import xsorb.ase_custom.xyz
from xsorb.ase_custom.atoms import AtomsCustom
from xsorb.ase_custom.xyz import write_xyz_custom

# Runtime patch for rendering is not called here but only explicitly
# with rom xsorb.ase_custom import povray when necessary
