import sys
from ase.io import read
from ase.optimize import BFGS

try:
    from ml_calculator_loader import NNloader
except:
    ImportError("ml_calculator_loader.py not found in the current folder.")

try:
    calculator = NNloader()
except:
    RuntimeError("Error loading the ML calculator.")


in_file = sys.argv[1]
traj_file = sys.argv[2]
log_file = sys.argv[3]
fmax = float(sys.argv[4])

structure = read(in_file)
structure.set_calculator(calculator)

# Optimize the structure
opt = BFGS(structure, logfile=log_file, trajectory=traj_file, maxstep=0.1)
opt.run(fmax=fmax, steps=1000) #limit the number of steps: if it doesn't converge in 1000 steps, it will never converge

with open(log_file, "a") as f:
    f.write("\nOptimization completed.\n")
