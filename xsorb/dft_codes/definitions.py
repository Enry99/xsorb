from ase.units import create_units
from xsorb.common_definitions import screening_outdir, relax_outdir, preopt_outdir

SUPPORTED_PROGRAMS = ['VASP', 'ESPRESSO']

HYBRID_SCREENING_THRESHOLDS = {
    'VASP' : -0.5, # ~ -2e-2 Ry/Bohr
    'ESPRESSO' : [5e-3, 5e-2]
}

UNITS_TO_EV_FACTOR = {
    'VASP' : 1,
    'ESPRESSO': create_units('2006')['Rydberg'],
    'ML': 1
}

OPTIMIZATION_COMPLETED_STRINGS = {
        #'VASP' : 'reached required accuracy - stopping structural energy minimisation', #in OUTCAR
        'VASP' : 'finalpos', #in vasprun.xml
        'ESPRESSO': 'Begin final coordinates',
        'ML': 'Optimization completed.'
}


SCF_NONCONVERGED_STRINGS = {
    'VASP': 'abcdefgxyz', #TODO: vasp does not stop the relax if one scf does not converge, so not necessary
    'ESPRESSO': 'convergence NOT achieved'
}

SCF_CONVERGED_STRINGS = {
    'VASP': 'abcdefgxyz', #TODO: update
    'ESPRESSO': '!'
}


OUT_FILE_PATHS = {
    'SCREENING': {
        'VASP': screening_outdir+'/{0}/vasprun.xml',
        'ESPRESSO': screening_outdir+'/{0}/screening_{0}.pwo',
    },

    'RELAX': {
        'VASP': relax_outdir+'/{0}/vasprun.xml',
        'ESPRESSO': relax_outdir+'/{0}/relax_{0}.pwo',
    },
    
    'PREOPT': {
        'ML': preopt_outdir+'/{0}/preopt_{0}.traj'
    }
}

IN_FILE_PATHS = {
    'SCREENING': {
        'VASP': screening_outdir+'/{0}/POSCAR',
        'ESPRESSO': screening_outdir+'/{0}/screening_{0}.pwi',
    },

    'RELAX': {
        'VASP': relax_outdir+'/{0}/POSCAR',
        'ESPRESSO': relax_outdir+'/{0}/relax_{0}.pwi',
    },

    'PREOPT': {
        'ML': preopt_outdir+'/{0}/preopt_{0}.xyz'
    }
}

LOG_FILE_PATHS = {
    'SCREENING': {
        'VASP': screening_outdir+'/{0}/vasprun.xml',
        'ESPRESSO': screening_outdir+'/{0}/screening_{0}.pwo',
    },

    'RELAX': {
        'VASP': relax_outdir+'/{0}/vasprun.xml',
        'ESPRESSO': relax_outdir+'/{0}/relax_{0}.pwo',
    },

    'PREOPT': {
        'ML' : preopt_outdir+'/{0}/preopt_{0}.log'
    }
}

SBATCH_POSTFIX = {
    'SCREENING': {
        'VASP': '',
        'ESPRESSO': 'screening_{0}.pwi screening_{0}.pwo',
    },

    'RELAX': {
        'VASP': '',
        'ESPRESSO': 'relax_{0}.pwi relax_{0}.pwo',
    },

}


FRAGMENTS_IN_FILE_PATHS = {
    'VASP': 'fragments/{0}/POSCAR',
    'ESPRESSO': 'fragments/{0}/{0}.pwi',
    'ML': 'fragments/{0}/{0}_ml.xyz'
}

FRAGMENTS_OUT_FILE_PATHS = {
    'VASP': 'fragments/{0}/vasprun.xml',
    'ESPRESSO': 'fragments/{0}/{0}.pwo',
    'ML': 'fragments/{0}/{0}_ml.traj'
}

FRAGMENTS_LOG_FILE_PATHS = {
    'VASP': 'fragments/{0}/vasprun.xml',
    'ESPRESSO': 'fragments/{0}/{0}.pwo',
    'ML': 'fragments/{0}/{0}_ml.log'
}

SBATCH_POSTFIX_FRAGS = {
    'VASP': '',
    'ESPRESSO': '{0}.pwi {0}.pwo',
}