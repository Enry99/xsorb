import sys, os
from ase.io import read
from ase.optimize import BFGSLineSearch
import numpy as np


in_file = sys.argv[1]
out_file = sys.argv[2]
log_file = sys.argv[3]
main_dir = sys.argv[4]
fix = sys.argv[5] == 'fix'
slab_ibeg = int(sys.argv[6])
slab_iend = int(sys.argv[7])

try:
    sys.path.append(os.path.abspath(main_dir))
    from ml_calculator_loader import NNloader
except:
    raise ImportError("ml_calculator_loader.py not found in the current folder.")

try:
    calculator = NNloader()
except:
    raise RuntimeError("Error loading the ML calculator.")


atoms = read(in_file)
atoms.calc = calculator

if fix:
    from ase.neighborlist import NeighborList, natural_cutoffs
    from ase.constraints import FixBondLengths
    cutoffs = natural_cutoffs(atoms, mult=1.1)

    slab_indices = [atom.index for atom in atoms if slab_ibeg <= atom.index <= slab_iend]
    mol_indices = [atom.index for atom in atoms if atom.index not in slab_indices]            
    nl = NeighborList(cutoffs, skin=0, self_interaction=False)
    nl.update(atoms)
    cm = nl.get_connectivity_matrix()

    slab_connected_pairs = [(i,j) for i in slab_indices for j in slab_indices if cm[i,j]]
    mol_connected_pairs = [(i,j) for i in mol_indices for j in mol_indices if cm[i,j]]

    atoms.set_constraint([FixBondLengths(slab_connected_pairs, tolerance=0.2), FixBondLengths(mol_connected_pairs, tolerance=0.3)])
    print("Constrained optimization.")


# Optimize the structure
opt = BFGSLineSearch(atoms, trajectory=out_file, logfile=log_file, maxstep=0.1)
opt.run(fmax=0.05, steps=500) #limit the number of steps: if it doesn't converge in 500 steps, it will never converge

with open(log_file, "a") as f:
    f.write("\nOptimization completed.\n")