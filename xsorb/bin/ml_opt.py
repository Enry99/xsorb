#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Enrico Pedretti

'''
Script to optimize a structure using a ML calculator.
'''

#TODO: when creating the package, make sure that this script
# is available as xsorb-ml-opt

import sys
import os

from ase.optimize import BFGS, BFGSLineSearch

from xsorb.io.utils import ase_custom_read as read

def main():
    '''
    Command-line script to optimize a structure using a ML calculator.
    Usage: python3 ml_opt.py in_file out_file log_file main_dir
    '''

    if len(sys.argv) != 5:
        raise ValueError("Usage: python3 in_file out_file log_file main_dir")

    in_file = sys.argv[1]
    out_file = sys.argv[2]
    log_file = sys.argv[3]
    main_dir = sys.argv[4] #path to the directory containing the ml_calculator_loader.py


    try:
        sys.path.append(os.path.abspath(main_dir))
        from ml_calculator_loader import NNloader
    except ImportError as exc:
        raise ImportError("ml_calculator_loader.py not found.") from exc

    try:
        calculator = NNloader()
    except Exception as exc:
        raise RuntimeError("Error loading the ML calculator.") from exc


    atoms = read(in_file)
    atoms.calc = calculator


    # Optimize the structure

    #First, try with BFGSLinesearch for 500 steps
    opt = BFGSLineSearch(atoms, trajectory=out_file, logfile=log_file, maxstep=0.1)
    converged = opt.run(fmax=0.01, steps=500)

    #If not converged, try with regular BFGS
    if not converged:
        opt = BFGS(atoms, trajectory=out_file, logfile=log_file, append_trajectory=True)
        converged = opt.run(fmax=0.01, steps=500)

    with open(log_file, "a", encoding=sys.getfilesystemencoding()) as f:
        f.write(f"\nOptimization {'converged' if converged else 'NOT converged'}.\n")


if __name__ == '__main__':
    main()
