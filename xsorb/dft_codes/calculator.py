#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Enrico Pedretti

"""
DFT CODE (and ML)-specific functions
to write the input files for the calculations

"""

from __future__ import annotations
from typing import Callable
import shutil
import glob
import os
import sys
import warnings
from pathlib import Path

from pymatgen.io.vasp.sets import (MPRelaxSet, MPMetalRelaxSet, MPScanRelaxSet,
                                   MPHSERelaxSet, MITRelaxSet)
from pymatgen.io.ase import AseAtomsAdaptor
from ase import Atoms
from ase.constraints import FixScaled
from ase.calculators.espresso import Espresso
from ase.calculators.vasp import Vasp

from xsorb.ase_custom.io import write


class MLFakeCalculator():
    '''
    Fake calculator for the ML optimization,
    which just writes the input files and does not perform any calculation.
    Its structure is made to be compatible with the other calculators,
    so that the write_input function can be called in the same way,
    and has the same initialization parameters.
    '''

    def __init__(self, label : str, directory : str):
        self.label = label
        self.directory = Path(directory)


    def write_input(self, atoms : Atoms):
        '''
        Function that writes the input file for the ML optimization
        Works in the same way as the write_input function of the other calculators,
        creating the directory if it does not exist
        '''
        self.directory.mkdir(exist_ok=True, parents=True)
        write(self.directory / f'{self.label}.xyz', atoms)


def setup_ML_calculator(**kwargs) -> MLFakeCalculator:
    '''
    Setup the ML calculator for the optimization and returns the calculator,
    to be used in a code-agnostic way in the write_file_with_calculator function.
    Must have this signature to be compatible with the other calculators.
    '''

    label : str = kwargs.get("label")
    directory : str = kwargs.get("directory")

    return MLFakeCalculator(label, directory)


def setup_Espresso_calculator(**kwargs) -> Espresso:
    '''
    Setup the Espresso calculator for the DFT calculation
    and returns the calculator, to be used in a code-agnostic way in
    the write_file_with_calculator function.
    '''

    dftsettings : dict = kwargs.get("dftsettings")
    label : str = kwargs.get("label")
    directory : str = kwargs.get("directory")

    return Espresso(label = label,
                    directory=directory,
                    pseudopotentials=dftsettings['pseudopotentials'],
                    kpts=dftsettings["kpts"],
                    koffset=dftsettings["koffset"],
                    input_data=dftsettings,
                    additional_cards=dftsettings['additional_cards'])


def setup_Vasp_calculator(**kwargs) -> Vasp:
    '''
    Setup the Vasp calculator for the DFT calculation
    and returns the calculator, to be used in a code-agnostic way in
    the write_file_with_calculator function.
    '''

    atoms : Atoms = kwargs.get("atoms")
    dftsettings : dict = kwargs.get("dftsettings")
    directory : str = kwargs.get("directory")

    preset_incar_settings = {}

    if "pymatgen_set" in dftsettings:
        sets_map = {'mprelaxset': MPRelaxSet,
                    'mpmetalrelaxset': MPMetalRelaxSet,
                    'mpscanrelaxset': MPScanRelaxSet,
                    'mphserelaxset': MPHSERelaxSet,
                    'mitrelaxset': MITRelaxSet}
        RelaxSet = sets_map[dftsettings["pymatgen_set"]]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore") #to suppress constraints not supported in pymatgen
            relax = RelaxSet(AseAtomsAdaptor.get_structure(atoms))

        preset_incar_settings = {k.lower(): v for k, v in relax.incar.as_dict().items()}
        preset_incar_settings.pop('@module')
        preset_incar_settings.pop('@class')

        if "kpoints_string" not in dftsettings:
            relax.kpoints.write_file('_temp_kpts_')

    preset_incar_settings["xc"] = dftsettings["vasp_xc_functional"]

    adjust_constraints(atoms, 'vasp')

    os.environ["VASP_PP_PATH"] = dftsettings["vasp_pp_path"]

    calc = Vasp(directory=directory,
                setups=dftsettings["vasp_pseudo_setups"],
                **preset_incar_settings) #set here the default values from pymat recommended


    #write user-defined settings to string, to be parsed by ASE, overriding the preset flags
    if "incar_string" in dftsettings:
        with open('_temp_incar_', 'w',encoding=sys.getfilesystemencoding()) as f:
            f.write(dftsettings["incar_string"])
        calc.read_incar('_temp_incar_')
        os.remove('_temp_incar_')

    if "kpoints_string" in dftsettings:
        with open('_temp_kpts_', 'w',encoding=sys.getfilesystemencoding()) as f:
            f.write(dftsettings["kpoints_string"])

    if "kpoints_string" in dftsettings or "pymatgen_set" in dftsettings:
        calc.read_kpoints('_temp_kpts_')
        os.remove('_temp_kpts_')

    return calc


def write_file_with_calculator(atoms : Atoms,
                               program: str,
                               dftsettings : dict,
                               label : str,
                               directory : str):
    '''
    Write the input files for the DFT calculation, using the specified calculator.
    '''

    setup_functions : dict[str, Callable[...,Vasp|Espresso|MLFakeCalculator]] = {
        'espresso': setup_Espresso_calculator,
        'vasp': setup_Vasp_calculator,
        'ml': setup_ML_calculator}


    if program not in setup_functions:
        raise ValueError(f'Program {program} not supported')


    calc = setup_functions[program](atoms, dftsettings, label, directory)
    calc.write_input(atoms)


def adjust_constraints(atoms : Atoms, program : str):
    '''
    Convert the constraints from the FixScaled format to the FixCartesian format
    for the VASP calculator. The constraints are set in the Atoms object (inplace)
    '''

    if program == 'vasp':
        #convert from FixCartesian to FixScaled
        c = [FixScaled(a=constr.a, mask=constr.mask) \
             for constr in atoms.constraints]
        atoms.set_constraint(c)


def edit_files_for_restart(program : str, paths : list[str]):
    '''
    Edit the input files, setting the correct flags for restart.

    Args:
    - program: DFT program. Possible values: 'espresso' or 'vasp'
    - paths: list of paths of the input files to be edited
    '''

    for path in paths:
        if program == 'espresso':
            with open(path, 'r',encoding=sys.getfilesystemencoding()) as f:
                lines = f.readlines()
                for i, line in enumerate(lines):
                    if 'from_scratch' in line:
                        lines[i] = lines[i].replace('from_scratch','restart')
                        break
            with open(path, 'w',encoding=sys.getfilesystemencoding()) as f:
                f.writelines(lines)

        elif program == 'vasp':
            poscar = path
            contcar = poscar.replace('POSCAR', 'CONTCAR')
            incar = poscar.replace('POSCAR', 'INCAR')
            outcar = poscar.replace('POSCAR', 'OUTCAR')
            vasprun = poscar.replace('POSCAR', 'vasprun.xml')
            oszicar = poscar.replace('POSCAR', 'OSZICAR')

            with open(incar, 'r',encoding=sys.getfilesystemencoding()) as f:
                lines = f.readlines()
                istart_found = False
                for i, line in enumerate(lines):
                    if 'ISTART' in line:
                        lines[i] = 'ISTART = 1\n'
                        istart_found = True
                        break
                if not istart_found: lines.append('ISTART = 1\n')

            with open(incar, 'w',encoding=sys.getfilesystemencoding()) as f:
                f.writelines(lines)

            #copy files for the first part of the relaxation, to avoid overwriting
            last_i = len( glob.glob( outcar.replace('OUTCAR', 'OUTCAR_') ) )
            shutil.copyfile(outcar, outcar.replace('OUTCAR', f'OUTCAR_{last_i}'))
            shutil.copyfile(vasprun, vasprun.replace('vasprun.xml', f'vasprun_{last_i}.xml'))
            shutil.copyfile(poscar, poscar.replace('POSCAR', f'POSCAR_{last_i}'))
            shutil.copyfile(oszicar, oszicar.replace('OSZICAR', f'OSZICAR_{last_i}'))

            #copy contcar to poscar to restart
            shutil.copyfile(contcar, poscar)
