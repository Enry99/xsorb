"""
@author: Enrico Pedretti

DFT CODE-specific functions (all others must be code-agnostic)

"""

import shutil, glob, os
from ase.constraints import FixCartesian, FixScaled
from ase.calculators.espresso import Espresso
from ase.calculators.vasp import Vasp
from settings import Settings
from filenames import *


UNITS_TO_EV_FACTOR = {
    'VASP' : 1,
    'ESPRESSO': rydbergtoev
}

COMPLETION_STRINGS = {
    'RELAX_COMPLETED': {
        #'VASP' : 'reached required accuracy - stopping structural energy minimisation', #in OUTCAR
        'VASP' : 'finalpos', #in vasprun.xml
        'ESPRESSO': 'Final energy'
    },

    'SCF_NONCONVERGED': {
        'VASP': 'abcdefgxyz', #TODO: vasp does not stop the relax if one scf does not converge, so not necessary
        'ESPRESSO': 'convergence NOT achieved'
    }
}

OUT_FILE_PATHS = {
    'SCREENING': {
        'VASP': screening_outdir+'/{0}/vasprun.xml',
        'ESPRESSO': 'screening_{0}.pwo',
    },

    'RELAX': {
        'VASP': relax_outdir+'/{0}/vasprun.xml',
        'ESPRESSO': 'relax_{0}.pwo',
    }
}

IN_FILE_PATHS = {
    'SCREENING': {
        'VASP': screening_outdir+'/{0}/POSCAR',
        'ESPRESSO': 'screening_{0}.pwi',
    },

    'RELAX': {
        'VASP': relax_outdir+'/{0}/POSCAR',
        'ESPRESSO': 'relax_{0}.pwi',
    }
}

SBATCH_POSTFIX = {
    'SCREENING': {
        'VASP': '',
        'ESPRESSO': '{0}/screening_{1}.pwi {0}/screening_{1}.pwo',
    },

    'RELAX': {
        'VASP': '',
        'ESPRESSO': '{0}/screening_{1}.pwi {0}/screening_{1}.pwo',
    }
}




def override_settings(settings : Settings, calc_type : str):
    '''
    Overrides the DFT code settings with specific options for the preliminary screening or the final relaxation

    Args:
    - settings: Settings object, containing the dft code parameters
    - calc_type: 'SCREENING' or 'RELAX'
    '''
    
    if calc_type is 'SCREENING':
        if settings.program is 'ESPRESSO':
            settings.espresso_settings_dict['CONTROL'].update({'calculation' : 'relax' })
            settings.espresso_settings_dict['CONTROL'].update({'etot_conv_thr' : settings.screening_conv_thr[0]})
            settings.espresso_settings_dict['CONTROL'].update({'forc_conv_thr' : settings.screening_conv_thr[1]})
            settings.espresso_settings_dict['IONS'].update({'upscale': 1})
            settings.espresso_settings_dict['CONTROL'].update({'restart_mode' : 'from_scratch'})

        elif settings.program is 'VASP':
            
            s = settings.incar_string.split('\n')
            
            missing_ibrion = True
            for i, line in enumerate(s):
                if 'EDIFFG' in line:
                    s[i] = f'EDIFFG = {settings.screening_conv_thr[0]}'
                if 'IBRION' in line: missing_ibrion = False
            if missing_ibrion: s.append('IBRION = 2')

            settings.incar_string = '\n'.join(s)


    elif calc_type is 'RELAX':
        if settings.program is 'ESPRESSO':
            settings.espresso_settings_dict['CONTROL'].update({'calculation' : 'relax'})
            settings.espresso_settings_dict['CONTROL'].update({'restart_mode' : 'from_scratch'})
            settings.espresso_settings_dict['IONS'].update({'ion_dynamics': settings.ion_dynamics})
        
        elif settings.program is 'VASP':
            s = settings.incar_string.split('\n')
            
            missing_ibrion = True
            for i, line in enumerate(s):
                if 'IBRION' in line: missing_ibrion = False
            if missing_ibrion: s.append('IBRION = 2')

            settings.incar_string = '\n'.join(s)


def Calculator(settings : Settings, label : str, atoms, directory : str):
    
    if settings.program is 'ESPRESSO':
        return Espresso(label = label, **settings.espresso_settings_dict)
    
    elif settings.program is 'VASP':
        
        preset_incar_settings = {}
        if settings.pymatgen_set:
            from pymatgen.io.vasp.sets import MPRelaxSet, MPMetalRelaxSet, MPScanRelaxSet, MPHSERelaxSet, MITRelaxSet
            from pymatgen.io.ase import AseAtomsAdaptor

            if settings.pymatgen_set.lower() is 'MPRelaxSet'.lower(): RelaxSet = MPRelaxSet
            elif settings.pymatgen_set.lower() is 'MPMetalRelaxSet'.lower(): RelaxSet = MPMetalRelaxSet
            elif settings.pymatgen_set.lower() is 'MPScanRelaxSet'.lower(): RelaxSet = MPScanRelaxSet
            elif settings.pymatgen_set.lower() is 'MPHSERelaxSet'.lower(): RelaxSet = MPHSERelaxSet            
            elif settings.pymatgen_set.lower() is 'MITRelaxSet'.lower(): RelaxSet = MITRelaxSet
            else: raise ValueError('Pymatgen preset not recognized.')

            atomscopy = atoms.copy()
            del atomscopy.constraints #to suppress the warning about constraints not supported in pymatgen
            relax = RelaxSet(AseAtomsAdaptor.get_structure(AseAtomsAdaptor.get_structure(atomscopy)))
            preset_incar_settings = {k.lower(): v for k, v in relax.incar.as_dict().items()} 
            preset_incar_settings.pop('@module')
            preset_incar_settings.pop('@class')
            relax.kpoints.write_file('_temp_kpts_')


        adjust_constraints(atoms, 'VASP')
        os.environ["VASP_PP_PATH"] = settings.vasp_pp_path
        calc = Vasp(directory=directory, 
                    xc=settings.vasp_xc_functional, 
                    setups=settings.vasp_pseudo_setups, 
                    **preset_incar_settings) #set here the default values from pymat recommended


        #write user-defined settings to string, to be parsed by ASE
        if settings.incar_string:
            with open('_temp_incar_', 'w') as f:
                f.write(settings.incar_string)

        if settings.kpoints_string:
            with open('_temp_kpts_', 'w') as f:
                f.write(settings.kpoints_string)

        calc.read_incar('_temp_incar_') 
        calc.read_kpoints('_temp_kpts_')
        
        return calc



def adjust_constraints(atoms, program : str):
    
    if program is 'VASP':
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
        if program is 'ESPRESSO':
            with open(IN_FILE_PATHS[calc_type][program].format(index), 'rw') as f:
                lines = f.readlines()
                for i, line in enumerate(lines):
                    if 'from_scratch' in line:
                        lines[i] = lines[i].replace('from_scratch','restart')
                        break
                f.writelines(lines)

        elif program is 'VASP':
            poscar = IN_FILE_PATHS[calc_type][program].format(index)
            contcar = poscar.replace('POSCAR', 'CONTCAR')
            incar = poscar.replace('POSCAR', 'INCAR')
            outcar = poscar.replace('POSCAR', 'OUTCAR')
            vasprun = poscar.replace('POSCAR', 'vasprun.xml')
            oszicar = poscar.replace('POSCAR', 'OSZICAR')
            
            with open(incar, 'rw') as f:
                lines = f.readlines()
                for i, line in enumerate(lines):
                    if 'ISTART' in line:
                        lines[i] = 'ISTART = 1'
                        break
                f.writelines(lines)

            #copy files for the first part of the relaxation, to avoid overwriting
            last_i = len( glob.glob( outcar.replace('OUTCAR', 'OUTCAR*') ) ) - 1
            shutil.copyfile(outcar, outcar.replace('OUTCAR', f'OUTCAR_{last_i}'))
            shutil.copyfile(vasprun, vasprun.replace('vasprun.xml', f'vasprun_{last_i}.xml'))
            shutil.copyfile(poscar, poscar.replace('POSCAR', f'POSCAR_{last_i}'))
            shutil.copyfile(oszicar, oszicar.replace('OSZICAR', f'OSZICAR_{last_i}'))
           
            #copy contcar to poscar to restart
            shutil.copyfile(contcar, poscar)




#from ase.calculators.vasp
def read_convergence(self, lines=None):
    """Method that checks whether a calculation has converged."""
    if not lines:
        lines = self.load_file('OUTCAR')

    converged = None
    # First check electronic convergence
    for line in lines:
        if 0:  # vasp always prints that!
            if line.rfind('aborting loop') > -1:  # scf failed
                raise RuntimeError(line.strip())
                break
        if 'EDIFF  ' in line:
            ediff = float(line.split()[2])
        if 'total energy-change' in line:
            # I saw this in an atomic oxygen calculation. it
            # breaks this code, so I am checking for it here.
            if 'MIXING' in line:
                continue
            split = line.split(':')
            a = float(split[1].split('(')[0])
            b = split[1].split('(')[1][0:-2]
            # sometimes this line looks like (second number wrong format!):
            # energy-change (2. order) :-0.2141803E-08  ( 0.2737684-111)
            # we are checking still the first number so
            # let's "fix" the format for the second one
            if 'e' not in b.lower():
                # replace last occurrence of - (assumed exponent) with -e
                bsplit = b.split('-')
                bsplit[-1] = 'e' + bsplit[-1]
                b = '-'.join(bsplit).replace('-e', 'e-')
            b = float(b)
            if [abs(a), abs(b)] < [ediff, ediff]:
                converged = True
            else:
                converged = False
                continue
    # Then if ibrion in [1,2,3] check whether ionic relaxation
    # condition been fulfilled
    if ((self.int_params['ibrion'] in [1, 2, 3]
            and self.int_params['nsw'] not in [0])):
        if not self.read_relaxed():
            converged = False
        else:
            converged = True
    return converged

def read_relaxed(self, lines=None):
    """Check if ionic relaxation completed"""
    if not lines:
        lines = self.load_file('OUTCAR')
    for line in lines:
        if 'reached required accuracy' in line:
            return True
    return False