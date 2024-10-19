#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Enrico Pedretti

'''
Constants related to DFT codes and their usage in the workflow,
such as file paths, completion checks, etc.
'''

SUPPORTED_PROGRAMS = ['vasp', 'espresso', 'ml']

SCREENING_OUTDIR            = 'screening_outdirs'
RELAX_OUTDIR                = 'relax_outdirs'
ML_OPT_OUTDIR               = 'ml_opt_outdirs'


HYBRID_SCREENING_THRESHOLDS = {
    'vasp' : -0.5,              # eV/A, ~ -2e-2 Ry/Bohr
    'espresso' : (5e-3, 5e-2)   # [Ry, Ry/Bohr]
}


# File paths #######################################################
#use pwi also for ml since they retain constraint, while xyz does not

IN_FILE_PATHS = {
    'screening': {
        'vasp': SCREENING_OUTDIR+'/{0}/POSCAR',
        'espresso': SCREENING_OUTDIR+'/{0}/screening_{0}.pwi',
    },

    'relax': {
        'vasp': RELAX_OUTDIR+'/{0}/POSCAR',
        'espresso': RELAX_OUTDIR+'/{0}/relax_{0}.pwi',
    },

    'preopt': {
        'ml': ML_OPT_OUTDIR+'/{0}/preopt_{0}.pwi'
    },

    'slab': {
        'vasp': 'slab/DFT/POSCAR',
        'espresso': 'slab/DFT/slab.pwi',
        'ml': 'slab/ML/slab.pwi'
    },

    'mol': {
        'vasp': 'mol/DFT/POSCAR',
        'espresso': 'mol/DFT/mol.pwi',
        'ml': 'mol/ML/mol.pwi'
    },

}

OUT_FILE_PATHS = {
    'screening': {
        'vasp': SCREENING_OUTDIR+'/{0}/vasprun.xml',
        'espresso': SCREENING_OUTDIR+'/{0}/screening_{0}.pwo',
    },

    'relax': {
        'vasp': RELAX_OUTDIR+'/{0}/vasprun.xml',
        'espresso': RELAX_OUTDIR+'/{0}/relax_{0}.pwo',
    },

    'preopt': {
        'ml': ML_OPT_OUTDIR+'/{0}/preopt_{0}.traj'
    },

    'slab': {
        'vasp': 'slab/DFT/vasprun.xml',
        'espresso': 'slab/DFT/slab.pwo',
        'ml': 'slab/ML/slab.traj'
    },

    'mol': {
        'vasp': 'mol/DFT/vasprun.xml',
        'espresso': 'mol/DFT/mol.pwo',
        'ml': 'mol/ML/mol.traj'
    }
}

LOG_FILE_PATHS = {
    'screening': {
        'vasp': SCREENING_OUTDIR+'/{0}/vasprun.xml',
        'espresso': SCREENING_OUTDIR+'/{0}/screening_{0}.pwo',
    },

    'relax': {
        'vasp': RELAX_OUTDIR+'/{0}/vasprun.xml',
        'espresso': RELAX_OUTDIR+'/{0}/relax_{0}.pwo',
    },

    'preopt': {
        'ml' : ML_OPT_OUTDIR+'/{0}/preopt_{0}.log'
    },

    'slab': {
        'vasp': 'slab/DFT/vasprun.xml',
        'espresso': 'slab/DFT/slab.pwo',
        'ml': 'slab/ML/slab.log'
    },

    'mol': {
        'vasp': 'mol/DFT/vasprun.xml',
        'espresso': 'mol/DFT/mol.pwo',
        'ml': 'mol/ML/mol.log'
    }
}

# Completion checks #########################################
OPTIMIZATION_COMPLETED_STRINGS = {
    #'vasp' : 'reached required accuracy - stopping structural energy minimisation', #in OUTCAR
    'vasp' : 'finalpos', #in vasprun.xml
    'espresso': 'Begin final coordinates',
    'ml': 'Optimization completed.'
}

SCF_NONCONVERGED_STRINGS = {
    'vasp': 'abcdefgxyz', #TODO: update
    'espresso': 'convergence NOT achieved'
}

SCF_CONVERGED_STRINGS = {
    'vasp': 'abcdefgxyz', #TODO: update
    'espresso': '!'
}


#Job submission #############################################
SBATCH_POSTFIX = {
    'vasp': '',
    'espresso': '{in_file} {out_file}',
    'ml': '{in_file} {out_file} {log_file} {main_dir}'
}

#TODO: update these
#Fragments ###################################################
FRAGMENTS_IN_FILE_PATHS = {
    'vasp': 'fragments/{0}/POSCAR',
    'espresso': 'fragments/{0}/{0}.pwi',
    'ml': 'fragments/{0}/{0}_ml.pwi'
}

FRAGMENTS_OUT_FILE_PATHS = {
    'vasp': 'fragments/{0}/vasprun.xml',
    'espresso': 'fragments/{0}/{0}.pwo',
    'ml': 'fragments/{0}/{0}_ml.traj'
}

FRAGMENTS_LOG_FILE_PATHS = {
    'vasp': 'fragments/{0}/vasprun.xml',
    'espresso': 'fragments/{0}/{0}.pwo',
    'ml': 'fragments/{0}/{0}_ml.log'
}

SBATCH_POSTFIX_FRAGS = {
    'vasp': '',
    'espresso': '{0}.pwi {0}.pwo',
}
