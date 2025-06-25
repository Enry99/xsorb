#Created on Fri 6 Feb 2023

"""
--------------------------------- Xsorb ------------------------------------------

Xsorb is a lightweight, flexible and user-friendly tool to generate and run
adsorption calculations of molecules on crystalline or amorphous surfaces.
The program can automatically generate initial adsorption configurations by
combining molecular rotations and adsorption sites, which can then be optimized
using DFT (with VASP or Quantum Espresso) or through machine learning potentials.

It employs a computationally efficient two(three)-step approach, with a first
*screening* step, in which all the initial adsorption configurations are optimized
with a larger convergence threshold, and a second *relax* step, where only a
selected subset of configurations is fully optimized.

Machine learning potentials can be employed:
i) to perform a fast pre-optmization before the ab initio screening step
ii) to replace the ab initio screening step, using DFT only for the final relax
iii) as a single-run tool to fastly explore the configuration space.

--------------------------------- Paper ------------------------------------------

If you use Xsorb in your research, pleas consider citing the following paper:
E. Pedretti, P. Restuccia, M.C. Righi, Comput. Phys. Commun. 291 (2023), 108827
https://doi.org/10.1016/j.cpc.2023.108827

------------------------------ Useful links --------------------------------------

Official repository:    https://gitlab.com/triboteam/xsorbed
Latest updates:         https://github.com/Enry99/xsorb
Documentation:          https://gitlab.com/triboteam/xsorbed/-/wikis/home

----------------------------------------------------------------------------------

"""

__version__ = "3.0beta"
