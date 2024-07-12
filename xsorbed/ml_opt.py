import sys, os
from ase.io import read
from ase.optimize import BFGS
from xsorbed.dftcode_specific import IN_FILE_PATHS, OUT_FILE_PATHS, LOG_FILE_PATHS


try:
    sys.path.append(os.getcwd())
    from ml_calculator_loader import NNloader
except:
    ImportError("ml_calculator_loader.py not found in the current folder.")

try:
    calculator = NNloader()
except:
    RuntimeError("Error loading the ML calculator.")


calc_idx = sys.argv[1]
in_file = IN_FILE_PATHS['PREOPT'].format(calc_idx)
out_file = OUT_FILE_PATHS['PREOPT'].format(calc_idx)
log_file = LOG_FILE_PATHS['PREOPT'].format(calc_idx)


atoms = read(in_file)
atoms.calc

# Optimize the structure
opt = BFGS(atoms, trajectory=out_file, logfile=log_file, maxstep=0.1)
opt.run(fmax=0.01, steps=1000) #limit the number of steps: if it doesn't converge in 1000 steps, it will never converge

with open(log_file, "a") as f:
    f.write("\nOptimization completed.\n")