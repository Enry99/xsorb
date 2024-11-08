#Created on Fri 6 Feb 2023

"""
--------------------------------- Xsorb ---------------------------------

Xsorb is a lightweight, flexible and user-friendly tool to generate and run ab initio calculations
of molecules adsorbed over substrates. The program can automatically generate initial adsorption
configurations combining molecular rotations and adsorption sites, and optimize them by ab initio
DFT calculations (with VASP or Quantum Espresso) or by a machine learning model.

It employs a computationally efficient two(three)-step approach, with a first *screening* step
in which all the initial adsorption configurations are optimized with a larger convergence threshold,
and a second *relax* step, where only a selected subset of configurations is fully optimized.

The machine learning model can be employed:
i) as a pre-optmization tool, using the final positions as an initial guess for the screening procedue
ii) to substitute the ab initio screening
iii) as a single-run tool to fastly explore the configuration space.

--------------------------------- Paper ---------------------------------

Please cite the following paper if you use Xsorb in your research:
E. Pedretti, P. Restuccia, M.C. Righi, Comput Phys Commun, 291 (2023), Article 108827
https://doi.org/10.1016/j.cpc.2023.108827

------------------------------ Useful links -----------------------------

Official repository:    https://gitlab.com/triboteam/xsorbed
Latest updates:         https://github.com/Enry99/xsorb
Documentation:          https://gitlab.com/triboteam/xsorbed/-/wikis/home

"""

__version__ = "3.0beta"
