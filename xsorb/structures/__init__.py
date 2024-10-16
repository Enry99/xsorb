"""
Contains the modules to generate the adsorption structures,
based on adsorption sites and molecule rotations

User-exposed classes and functions:
- AdsorptionStructuresGenerator: class to generate the adsorption structures
- slab_mol_bonds: function to get the bonds between a slab and a molecule

"""

from xsorb.structures.generation import AdsorptionStructuresGenerator, AdsorptionStructure
from xsorb.structures.slab import Slab
from xsorb.structures.molecule import Molecule
from xsorb.structures.utils import slab_mol_bonds, set_fixed_slab_constraints
