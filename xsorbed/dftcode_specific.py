"""
@author: Enrico Pedretti

DFT CODE-specific functions (all others must be code-agnostic)

"""

import shutil, glob
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
            pass

    elif calc_type is 'RELAX':
        pass


def Calculator(settings : Settings, label : str, atoms):
    if settings.program is 'ESPRESSO':
        return Espresso(label = label, **settings.espresso_settings_dict)
    if settings.program is 'VASP': 
        return Vasp()


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