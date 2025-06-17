'''
Module for defining file names and paths used in the xsorb package,
except for those related to dft codes, which are defined in xsorb.dft_codes.definitions
'''

STRUCTURES_DB_NAME = 'structures.json'

CALC_DB_NAMES = { # DO NOT CHANGE THE ORDER, as write_csvfile relies on it
    'mlopt': 'mlopt.json',
    'screening': 'screening.json',
    'relax': 'relaxations.json',
}

ALL_DB_NAMES = list(CALC_DB_NAMES.values()) + [STRUCTURES_DB_NAME]

ADSITES_FILENAME = 'adsites.json'

JOBS_FILENAME = '.submitted_jobs.txt'

SCREENING_OUTDIR            = 'screening_outdirs'
RELAX_OUTDIR                = 'relax_outdirs'
ML_OPT_OUTDIR               = 'ml_opt_outdirs'
SLAB_OUTDIR                = 'slab_outdirs'
MOL_OUTDIR                 = 'mol_outdirs'

# for convenience in cleanup
outdirs = [SCREENING_OUTDIR, RELAX_OUTDIR, ML_OPT_OUTDIR, SLAB_OUTDIR, MOL_OUTDIR]

# for convenience in check:
all_files_and_dirs = ALL_DB_NAMES + [ADSITES_FILENAME] + outdirs
