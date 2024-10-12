from ase.units import create_units
from xsorb.common_definitions import screening_outdir, relax_outdir, preopt_outdir

SUPPORTED_PROGRAMS = ['vasp', 'espresso']

HYBRID_SCREENING_THRESHOLDS = {
    'vasp' : -0.5, # ~ -2e-2 Ry/Bohr
    'espresso' : [5e-3, 5e-2]
}

UNITS_TO_EV_FACTOR = {
    'vasp' : 1,
    'espresso': create_units('2006')['Rydberg'],
    'ML': 1
}

OPTIMIZATION_COMPLETED_STRINGS = {
        #'vasp' : 'reached required accuracy - stopping structural energy minimisation', #in OUTCAR
        'vasp' : 'finalpos', #in vasprun.xml
        'espresso': 'Begin final coordinates',
        'ML': 'Optimization completed.'
}


SCF_NONCONVERGED_STRINGS = {
    'vasp': 'abcdefgxyz', #TODO: vasp does not stop the relax if one scf does not converge, so not necessary
    'espresso': 'convergence NOT achieved'
}

SCF_CONVERGED_STRINGS = {
    'vasp': 'abcdefgxyz', #TODO: update
    'espresso': '!'
}


OUT_FILE_PATHS = {
    'screening': {
        'vasp': screening_outdir+'/{0}/vasprun.xml',
        'espresso': screening_outdir+'/{0}/screening_{0}.pwo',
    },

    'relax': {
        'vasp': relax_outdir+'/{0}/vasprun.xml',
        'espresso': relax_outdir+'/{0}/relax_{0}.pwo',
    },

    'preopt': {
        'ML': preopt_outdir+'/{0}/preopt_{0}.traj'
    }
}

IN_FILE_PATHS = {
    'screening': {
        'vasp': screening_outdir+'/{0}/POSCAR',
        'espresso': screening_outdir+'/{0}/screening_{0}.pwi',
    },

    'relax': {
        'vasp': relax_outdir+'/{0}/POSCAR',
        'espresso': relax_outdir+'/{0}/relax_{0}.pwi',
    },

    'preopt': {
        'ML': preopt_outdir+'/{0}/preopt_{0}.xyz'
    }
}

LOG_FILE_PATHS = {
    'screening': {
        'vasp': screening_outdir+'/{0}/vasprun.xml',
        'espresso': screening_outdir+'/{0}/screening_{0}.pwo',
    },

    'relax': {
        'vasp': relax_outdir+'/{0}/vasprun.xml',
        'espresso': relax_outdir+'/{0}/relax_{0}.pwo',
    },

    'preopt': {
        'ML' : preopt_outdir+'/{0}/preopt_{0}.log'
    }
}

SBATCH_POSTFIX = {
    'screening': {
        'vasp': '',
        'espresso': 'screening_{0}.pwi screening_{0}.pwo',
    },

    'relax': {
        'vasp': '',
        'espresso': 'relax_{0}.pwi relax_{0}.pwo',
    },

}


FRAGMENTS_IN_FILE_PATHS = {
    'vasp': 'fragments/{0}/POSCAR',
    'espresso': 'fragments/{0}/{0}.pwi',
    'ML': 'fragments/{0}/{0}_ml.xyz'
}

FRAGMENTS_OUT_FILE_PATHS = {
    'vasp': 'fragments/{0}/vasprun.xml',
    'espresso': 'fragments/{0}/{0}.pwo',
    'ML': 'fragments/{0}/{0}_ml.traj'
}

FRAGMENTS_LOG_FILE_PATHS = {
    'vasp': 'fragments/{0}/vasprun.xml',
    'espresso': 'fragments/{0}/{0}.pwo',
    'ML': 'fragments/{0}/{0}_ml.log'
}

SBATCH_POSTFIX_FRAGS = {
    'vasp': '',
    'espresso': '{0}.pwi {0}.pwo',
}