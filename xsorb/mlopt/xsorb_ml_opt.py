#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Enrico Pedretti

'''
Script to optimize a structure using a ML calculator.
'''

import sys
import os

from ase.optimize import BFGSLineSearch, FIRE2

from xsorb.ase_custom.io import ase_custom_read as read

def main():
    '''
    Command-line script to optimize a structure using a ML calculator.
    Usage: python3 xsorb-ml-opt in_file out_file log_file main_dir
    '''

    if len(sys.argv) != 5:
        raise ValueError("Usage: xsorb-ml-opt in_file out_file log_file main_dir")

    in_file = sys.argv[1]
    out_file = sys.argv[2]
    log_file = sys.argv[3]
    main_dir = sys.argv[4] #path to the directory containing the ml_calculator_loader.py


    try:
        sys.path.append(os.path.abspath(main_dir))
        from ml_calculator_loader import NNloader #runtime import
    except ImportError as exc:
        raise ImportError("ml_calculator_loader.py not found.") from exc

    try:
        calculator = NNloader()
    except Exception as exc:
        raise RuntimeError("Error loading the ML calculator.") from exc


    atoms = read(in_file)
    atoms.calc = calculator


    # Optimize the structure

    #First, try with BFGSLinesearch for 300 steps
    opt = BFGSLineSearch(atoms, trajectory=out_file, logfile=log_file, maxstep=0.1)
    converged = opt.run(fmax=0.01, steps=300)

    #If not converged, try with FIRE (more robust, but can require more force calls)
    if not converged:
        opt = FIRE2(atoms, trajectory=out_file, logfile=log_file, append_trajectory=True)
        converged = opt.run(fmax=0.01, steps=2000)

    with open(log_file, "a", encoding=sys.getfilesystemencoding()) as f:
        f.write(f"\nOptimization {'converged' if converged else 'NOT converged'}.\n")
