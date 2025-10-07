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


HYBRID_SCREENING_THRESHOLD = 0.5  # force threshold (eV/,  ~2e-2 Ry/Bohr)
RELAX_THRESHOLD = 0.01            # force threshold (eV/A, ~4e-4 Ry/Bohr)


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

    'mlopt': {
        'ml': ML_OPT_OUTDIR+'/{0}/mlopt_{0}.xyz'
    },

    'slab': {
        'vasp': ISOLATED_OUTDIRS+'/DFT/slab/POSCAR',
        'espresso': ISOLATED_OUTDIRS+'/DFT/slab/slab.pwi',
        'ml': ISOLATED_OUTDIRS+'/ML/slab/slab.xyz'
    },

    'mol': {
        'vasp': ISOLATED_OUTDIRS+'/DFT/mol/POSCAR',
        'espresso': ISOLATED_OUTDIRS+'/DFT/mol/mol.pwi',
        'ml': ISOLATED_OUTDIRS+'/ML/mol/mol.xyz'
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

    'mlopt': {
        'ml': ML_OPT_OUTDIR+'/{0}/mlopt_{0}.traj'
    },

    'slab': {
        'vasp': ISOLATED_OUTDIRS+'/DFT/slab/vasprun.xml',
        'espresso': ISOLATED_OUTDIRS+'/DFT/slab/slab.pwo',
        'ml': ISOLATED_OUTDIRS+'/ML/slab/slab.traj'
    },

    'mol': {
        'vasp': ISOLATED_OUTDIRS+'/DFT/mol/vasprun.xml',
        'espresso': ISOLATED_OUTDIRS+'/DFT/mol/mol.pwo',
        'ml': ISOLATED_OUTDIRS+'/ML/mol/mol.traj'
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

    'mlopt': {
        'ml' : ML_OPT_OUTDIR+'/{0}/mlopt_{0}.log'
    },

    'slab': {
        'vasp': ISOLATED_OUTDIRS+'/DFT/slab/vasprun.xml',
        'espresso': ISOLATED_OUTDIRS+'/DFT/slab/slab.pwo',
        'ml': ISOLATED_OUTDIRS+'/ML/slab/slab.log'
    },

    'mol': {
        'vasp': ISOLATED_OUTDIRS+'/DFT/mol/vasprun.xml',
        'espresso': ISOLATED_OUTDIRS+'/DFT/mol/mol.pwo',
        'ml': ISOLATED_OUTDIRS+'/ML/mol/mol.log'
    }
}

# Completion checks #########################################
OPTIMIZATION_COMPLETED_STRINGS = {
    #'vasp' : 'reached required accuracy - stopping structural energy minimisation', #in OUTCAR
    'vasp' : 'finalpos', #in vasprun.xml
    'espresso': 'Begin final coordinates',
    'ml': 'Optimization converged.'
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
