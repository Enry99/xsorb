#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Enrico Pedretti

'''
Constants related to DFT codes and their usage in the workflow,
such as file paths, completion checks, etc.
'''

from xsorb.io.filenames import (
    SCREENING_OUTDIR,
    RELAX_OUTDIR,
    ML_OPT_OUTDIR,
    ISOLATED_OUTDIRS
)


SUPPORTED_PROGRAMS = ['vasp', 'espresso', 'ml']


HYBRID_SCREENING_THRESHOLDS = {
    'vasp' : -0.5,              # eV/A, ~ -2e-2 Ry/Bohr
    'espresso' : (5e-3, 5e-2)   # [Ry, Ry/Bohr]
}


# File paths #######################################################

IN_FILE_PATHS = {
    'screening': {
        'vasp': SCREENING_OUTDIR+'/{0}/POSCAR',
        'espresso': SCREENING_OUTDIR+'/{0}/screening_{0}.pwi',
    },

    'relax': {
        'vasp': RELAX_OUTDIR+'/{0}/POSCAR',
        'espresso': RELAX_OUTDIR+'/{0}/relax_{0}.pwi',
    },

    'mlopt': {
        'unified': ML_OPT_OUTDIR+'/{0}/mlopt_{0}.xyz'
    },

    'slab': {
        'vasp': ISOLATED_OUTDIRS+'/DFT/slab/POSCAR',
        'espresso': ISOLATED_OUTDIRS+'/DFT/slab/slab.pwi',
        'unified': ISOLATED_OUTDIRS+'/ML/slab/slab.xyz'
    },

    'mol': {
        'vasp': ISOLATED_OUTDIRS+'/DFT/mol/{0}/POSCAR',
        'espresso': ISOLATED_OUTDIRS+'/DFT/mol/{0}/mol.pwi',
        'unified': ISOLATED_OUTDIRS+'/ML/mol/{0}/mol.xyz'
    },

}

OUT_FILE_PATHS = {
    'screening': {
        'vasp': SCREENING_OUTDIR+'/{0}/vasprun.xml',
        'espresso': SCREENING_OUTDIR+'/{0}/screening_{0}.pwo',
        'unified': SCREENING_OUTDIR+'/{0}/screening_{0}.traj',
    },

    'relax': {
        'vasp': RELAX_OUTDIR+'/{0}/vasprun.xml',
        'espresso': RELAX_OUTDIR+'/{0}/relax_{0}.pwo',
        'unified': RELAX_OUTDIR+'/{0}/relax_{0}.traj',
    },

    'mlopt': {
        'unified': ML_OPT_OUTDIR+'/{0}/mlopt_{0}.traj'
    },

    'slab': {
        'vasp': ISOLATED_OUTDIRS+'/DFT/slab/vasprun.xml',
        'espresso': ISOLATED_OUTDIRS+'/DFT/slab/slab.pwo',
        'unified': ISOLATED_OUTDIRS+'/DFT/slab/slab.traj',
        'unified': ISOLATED_OUTDIRS+'/ML/slab/slab.traj'
    },

    'mol': {
        'vasp': ISOLATED_OUTDIRS+'/DFT/mol/{0}/vasprun.xml',
        'espresso': ISOLATED_OUTDIRS+'/DFT/mol/{0}/mol.pwo',
        'unified': ISOLATED_OUTDIRS+'/DFT/mol/{0}/mol.traj',
        'unified': ISOLATED_OUTDIRS+'/ML/mol/{0}/mol.traj'
    }
}

LOG_FILE_PATHS = {
    'screening': {
        'vasp': SCREENING_OUTDIR+'/{0}/vasprun.xml',
        'espresso': SCREENING_OUTDIR+'/{0}/screening_{0}.pwo',
        'unified': SCREENING_OUTDIR+'/{0}/screening_{0}.log',
    },

    'relax': {
        'vasp': RELAX_OUTDIR+'/{0}/vasprun.xml',
        'espresso': RELAX_OUTDIR+'/{0}/relax_{0}.pwo',
        'unified': RELAX_OUTDIR+'/{0}/relax_{0}.log',
    },

    'mlopt': {
        'unified': ML_OPT_OUTDIR+'/{0}/mlopt_{0}.log'
    },

    'slab': {
        'vasp': ISOLATED_OUTDIRS+'/DFT/slab/vasprun.xml',
        'espresso': ISOLATED_OUTDIRS+'/DFT/slab/slab.pwo',
        'unified': ISOLATED_OUTDIRS+'/DFT/slab/slab.log',
        'unified': ISOLATED_OUTDIRS+'/ML/slab/slab.log'
    },

    'mol': {
        'vasp': ISOLATED_OUTDIRS+'/DFT/mol/{0}/vasprun.xml',
        'espresso': ISOLATED_OUTDIRS+'/DFT/mol/{0}/mol.pwo',
        'unified': ISOLATED_OUTDIRS+'/DFT/mol/{0}/mol.log',
        'unified': ISOLATED_OUTDIRS+'/ML/mol/{0}/mol.log'
    }
}

# Completion checks #########################################
# OPTIMIZATION_COMPLETED_STRINGS = {
#     #'vasp' : 'reached required accuracy - stopping structural energy minimisation', #in OUTCAR
#     'vasp' : 'finalpos', #in vasprun.xml
#     'espresso': 'Begin final coordinates',
#     'unified': 'Optimization converged.',
#     'ml': 'Optimization converged.'
# }
OPTIMIZATION_COMPLETED_STRINGS = {
    'xml' : 'finalpos', #in vasprun.xml
    'pwo': 'Begin final coordinates', # in espresso pwo
    'log': 'Optimization converged.', # in ase log
}


SCF_NONCONVERGED_STRINGS = {
    'vasp': 'abcdefgxyz', #TODO: update
    'espresso': 'convergence NOT achieved'
}

SCF_CONVERGED_STRINGS = {
    'vasp': '', #TODO: update
    'espresso': '!'
}


#Job submission #############################################
RUN_LINE_POSTFIX = {
    'vasp': '',
    'espresso': '-in {in_file} >> {out_file}',
    'unified': ''
}

#TODO: update these
#Fragments ###################################################
FRAGMENTS_IN_FILE_PATHS = {
    'vasp': 'fragments/{0}/POSCAR',
    'espresso': 'fragments/{0}/{0}.pwi',
    'ml': 'fragments/{0}/{0}_ml.xyz'
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
