'''
Module for defining file names and paths used in the xsorb package,
except for those related to dft codes, which are defined in xsorb.dft_codes.definitions
'''

DB_NAMES = {
    'structures': 'structures.json',
    'mlopt': 'mlopt.json',
    'screening': 'screening.json',
    'relax': 'relaxations.json',
}

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
all_files_and_dirs = list(DB_NAMES.values()) + [ADSITES_FILENAME] + outdirs
