#Created on Fri 6 Feb 2023

"""
Xsorb is a lightweight, flexible and user-friendly tool to generate
and run ab initio calculations of molecules adsorbed over substrates,
relying on VASP and Quantum Espresso for DFT calculations and on
ASE & Pymatgen to generate the adsorption structures.

The program can automatically generate initial adsorption configurations combining
molecular rotations and adsorption sites, and optimize them by ab initio DFT calculations
(with VASP or Quantum Espresso) or by a machine learning model. It employs a
computationally efficient two(three)-step approach, with a first screening step in which
all the initial adsorption configurations are optimized with a larger convergence threshold,
and a second step in which only the most promising configurations are optimized with the
desired convergence threshold. The machine learning model can be used as a pre-optmization tool,
using the final positions as an initial guess for the screening procedue, or to substitute the
ab initio screening (or even as the only optimization tool).

Please cite the following paper if you use Xsorb in your research:
E. Pedretti, P. Restuccia, M.C. Righi, Comput Phys Commun, 291 (2023), Article 108827,
https://doi.org/10.1016/j.cpc.2023.108827

The official repository is hosted at:
https://gitlab.com/triboteam/xsorbed

and at the developer's personal repository (for the latest updates):
https://github.com/Enry99/xsorb

The documentation can be found at:
https://gitlab.com/triboteam/xsorbed/-/wikis/home

"""

__version__ = "2.2"
