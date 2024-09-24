"""
@author: Enrico Pedretti

DFT CODE-specific functions (all others must be code-agnostic)

"""

import shutil
import glob
import os
import warnings
from ase import Atoms
from ase.constraints import FixScaled
from ase.calculators.espresso import Espresso
from ase.calculators.vasp import Vasp
from xsorb.settings import Settings
from xsorb.dft_codes.definitions import IN_FILE_PATHS

#TODO: check that when reading magmoms from input, the order is then changed correctly after resorting the poscar during write_inputs


def Calculator(settings : Settings, label : str, atoms : Atoms, directory : str):
    
    if settings.program == 'ESPRESSO':
        return Espresso(label = label,
                        pseudopotentials=settings.dftprogram_settings_dict['pseudopotentials'],
                        kpts=settings.dftprogram_settings_dict["kpts"],
                        koffset=settings.dftprogram_settings_dict["koffset"],
                        input_data=settings.dftprogram_settings_dict,
                        additional_cards=settings.dftprogram_settings_dict['additional_cards'])
                            
    elif settings.program == 'VASP':
        
        preset_incar_settings = {}
        
        if "pymatgen_set" in settings.dftprogram_settings_dict:
            from pymatgen.io.vasp.sets import MPRelaxSet, MPMetalRelaxSet, MPScanRelaxSet, MPHSERelaxSet, MITRelaxSet
            from pymatgen.io.ase import AseAtomsAdaptor
            
            sets_map = {'mprelaxset': MPRelaxSet, 'mpmetalrelaxset': MPMetalRelaxSet, 'mpscanrelaxset': MPScanRelaxSet, 'mphserelaxset': MPHSERelaxSet, 'mitrelaxset': MITRelaxSet}
            RelaxSet = sets_map[settings.dftprogram_settings_dict["pymatgen_set"]]
            with warnings.catch_warnings(): 
                warnings.simplefilter("ignore") #to suppress the warning about constraints not supported in pymatgen
                relax = RelaxSet(AseAtomsAdaptor.get_structure(atoms))

            preset_incar_settings = {k.lower(): v for k, v in relax.incar.as_dict().items()} 
            if 'isif' in preset_incar_settings: preset_incar_settings['isif'] = 2 #OVERRIDE: we do not want a vc-relax, unless explicitly specified in &INCAR
            preset_incar_settings.pop('@module')
            preset_incar_settings.pop('@class')
            relax.kpoints.write_file('_temp_kpts_')

        preset_incar_settings["xc"] = settings.dftprogram_settings_dict["vasp_xc_functional"]

        adjust_constraints(atoms, 'VASP')

        os.environ["VASP_PP_PATH"] = settings.dftprogram_settings_dict["vasp_pp_path"]

        calc = Vasp(directory=directory, 
                    setups=settings.dftprogram_settings_dict["vasp_pseudo_setups"], 
                    **preset_incar_settings) #set here the default values from pymat recommended


        #write user-defined settings to string, to be parsed by ASE, overriding the preset flags
        if "incar_string" in settings.dftprogram_settings_dict:
            with open('_temp_incar_', 'w') as f:
                f.write(settings.dftprogram_settings_dict["incar_string"])
            calc.read_incar('_temp_incar_')
            os.remove('_temp_incar_') 

        if "kpoints_string" in settings.dftprogram_settings_dict:
            with open('_temp_kpts_', 'w') as f:
                f.write(settings.dftprogram_settings_dict["kpoints_string"])

        if "kpoints_string" in settings.dftprogram_settings_dict or "pymatgen_set" in settings.dftprogram_settings_dict:
            calc.read_kpoints('_temp_kpts_') 
            os.remove('_temp_kpts_')
        
        return calc



def adjust_constraints(atoms : Atoms, program : str):
    
    if program == 'VASP':
        c = [FixScaled(atoms.cell, constr.a, ~constr.mask) for constr in atoms.constraints] #atoms.constraints are FixCartesian
        atoms.set_constraint(c)


def edit_files_for_restart(program : str, calc_type : str, indices : list):
    '''
    Edit the input files, setting the correct flags for restart.

    Args:
    - program: DFT program. Possible values: 'ESPRESSO' or 'VASP'
    - calc_type: 'SCREENING' or 'RELAX'
    - indices_list: indices of the calculations
    '''

    for index in indices:
        if program == 'ESPRESSO':
            with open(IN_FILE_PATHS[calc_type][program].format(index), 'r') as f:
                lines = f.readlines()
                for i, line in enumerate(lines):
                    if 'from_scratch' in line:
                        lines[i] = lines[i].replace('from_scratch','restart')
                        break
            with open(IN_FILE_PATHS[calc_type][program].format(index), 'w') as f:
                f.writelines(lines)

        elif program == 'VASP':
            poscar = IN_FILE_PATHS[calc_type][program].format(index)
            contcar = poscar.replace('POSCAR', 'CONTCAR')
            incar = poscar.replace('POSCAR', 'INCAR')
            outcar = poscar.replace('POSCAR', 'OUTCAR')
            vasprun = poscar.replace('POSCAR', 'vasprun.xml')
            oszicar = poscar.replace('POSCAR', 'OSZICAR')
            
            with open(incar, 'r') as f:
                lines = f.readlines()
                istart_found = False
                for i, line in enumerate(lines):
                    if 'ISTART' in line:
                        lines[i] = 'ISTART = 1\n'
                        istart_found = True
                        break
                if not istart_found: lines.append('ISTART = 1\n')

            with open(incar, 'w') as f:
                f.writelines(lines)

            #copy files for the first part of the relaxation, to avoid overwriting
            last_i = len( glob.glob( outcar.replace('OUTCAR', 'OUTCAR*') ) ) - 1
            shutil.copyfile(outcar, outcar.replace('OUTCAR', f'OUTCAR_{last_i}'))
            shutil.copyfile(vasprun, vasprun.replace('vasprun.xml', f'vasprun_{last_i}.xml'))
            shutil.copyfile(poscar, poscar.replace('POSCAR', f'POSCAR_{last_i}'))
            shutil.copyfile(oszicar, oszicar.replace('OSZICAR', f'OSZICAR_{last_i}'))
           
            #copy contcar to poscar to restart
            shutil.copyfile(contcar, poscar)
