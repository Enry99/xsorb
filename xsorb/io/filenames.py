'''
Module for defining file names and paths used in the xsorb package,
except for those related to dft codes/ml, which are defined in xsorb.dft_codes.definitions
'''

STRUCTURES_DB_NAME = 'structures.json'

CALC_DB_NAMES = { # DO NOT CHANGE THE ORDER, as write_csvfile relies on it
    'mlopt': 'mlopt.json',
    'screening': 'screening.json',
    'relax': 'relaxations.json',
}

JOBS_FILENAME = '.submitted_jobs.txt' # for slab/mol, as their job ids are not written in the db

SCREENING_OUTDIR            = 'screening_outdirs'
RELAX_OUTDIR                = 'relax_outdirs'
ML_OPT_OUTDIR               = 'mlopt_outdirs'
ISOLATED_OUTDIRS            = 'isolated_outdirs'

CONFORMERS_FILENAME         = 'molecule_conformers.xyz'


######## for convenience in check and cleanup ########
ALL_DB_NAMES = list(CALC_DB_NAMES.values()) + [STRUCTURES_DB_NAME]
ALL_OUTDIRS = [SCREENING_OUTDIR, RELAX_OUTDIR, ML_OPT_OUTDIR] #, ISOLATED_OUTDIRS]
ALL_FILES_AND_DIRS = ALL_DB_NAMES + ALL_OUTDIRS
