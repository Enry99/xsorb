'''
Module for geometry optimization through unified ASE interface with calculator
'''

import sys, os
import time
import json
from dataclasses import dataclass

from dacite import from_dict

from ase.calculators.socketio import SocketIOCalculator
from ase import Atoms

from xsorb.ase_custom.io import ase_custom_read as read
from xsorb.dft_codes.interactive_calculators import setup_calculator
from xsorb.io.filenames import UNIFIED_OPTIMIZATION_PARAMS


@dataclass
class OptimizationParameters: # pylint: disable=missing-class-docstring, too-many-instance-attributes
    code: str
    optimizer : str
    in_file : str
    out_file: str
    log_file: str
    command: str|None
    pseudo_dir: str|None
    main_dir: str
    calc_id: int
    fmax: float
    maxsteps: int = 500


def main():

    with open(UNIFIED_OPTIMIZATION_PARAMS) as f:
        parameters = json.load(f)

    opt_params = from_dict(data_class=OptimizationParameters, data=parameters)


    atoms : Atoms = read(opt_params.in_file)
    calc = setup_calculator(code=opt_params.code,
                            input_filename=opt_params.in_file,
                            command=opt_params.command,
                            pseudo_dir=opt_params.pseudo_dir,
                            main_dir=opt_params.main_dir,
                            calc_id=opt_params.calc_id,
                            atoms=atoms)
    atoms.calc = calc


    t0 = time.time()

    if opt_params.optimizer == 'bfgs_linesearch':
        from ase.optimize import BFGSLineSearch
        opt = BFGSLineSearch(atoms, trajectory=opt_params.out_file, logfile=opt_params.log_file,
                            c1=0.1,
                            restart='restart.json',
                            append_trajectory=True)
    elif opt_params.optimizer == 'quasi_newton':
        from ase.optimize import GoodOldQuasiNewton
        opt = GoodOldQuasiNewton(atoms, trajectory=opt_params.out_file, logfile=opt_params.log_file,
                                restart='restart.json',
                                append_trajectory=True)
    elif opt_params.optimizer == 'gpmin' or opt_params.optimizer == 'gpmin_ml':
        from ase.optimize import GPMin

        if opt_params.optimizer == 'gpmin_ml':
            from ase.optimize.gpmin.prior import CalculatorPrior
            try:
                sys.path.append(os.path.abspath(opt_params.main_dir))
                from gpmin_calculator_loader import NNloader #runtime import
            except ImportError as exc:
                raise ImportError("gpmin_calculator_loader.py not found.") from exc
            try:
                prior_calc = NNloader()
            except Exception as exc:
                raise RuntimeError("Error loading the ML calculator.") from exc

            prior = CalculatorPrior(atoms, prior_calc)
        else:
            prior = None
        opt = GPMin(atoms, trajectory=opt_params.out_file, logfile=opt_params.log_file,
                    prior=prior, update_hyperparams = True,
                    restart='restart.json',
                    append_trajectory=True)
    else:
        raise ValueError(f"Optimizer {opt_params.optimizer} not recognized. "
                         "Available options: bfgs_linesearch, quasi_newton, gpmin.")

    converged = opt.run(fmax=opt_params.fmax, steps=opt_params.maxsteps)

    #If not converged, try with FIRE (more robust, but can require more force calls)
    if not converged:
        from ase.optimize import FIRE2
        opt = FIRE2(atoms, trajectory=opt_params.out_file, logfile=opt_params.log_file,
                    restart='restart.json',
                    append_trajectory=True)
        converged = opt.run(fmax=opt_params.fmax, steps=opt_params.maxsteps)
    tf = time.time()


    with open(opt_params.log_file, "a", encoding=sys.getfilesystemencoding()) as f:
        f.write(f"\nFinal energy: {atoms.get_potential_energy():.3f} eV\n")
        f.write(f"n. of steps: {opt.nsteps}, elapsed time {tf - t0} s ({(tf - t0)/3600} h)\n")
        f.write(f"Optimization {'converged' if converged else 'NOT converged'}.\n")


    if isinstance(calc, SocketIOCalculator):
        calc.close()
