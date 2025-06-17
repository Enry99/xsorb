'''
Moudule to clean all files for an Xsorb run.
'''

import os
import shutil

from xsorb.io.filenames import DB_NAMES, ADSITES_FILENAME, JOBS_FILENAME, outdirs

def cleanup_xsorb_run(calc_only: bool = False):
    '''
    Clean all files for an Xsorb run.
    '''

    # ask the user for confirmation
    confirm = input("Do you really want to remove all files for this Xsorb run? (yes/no): ")
    if confirm != 'yes':
        print("Cleanup aborted.")
        return

    if calc_only:
        print("Cleaning up only calculation files...")

        for directory in outdirs:
            if os.path.exists(directory):
                shutil.rmtree(directory, ignore_errors=True)
                print(f"Removed directory: {directory}")
        print("Calculation files cleaned up.")

        return

    # Remove all database files
    for db_name in DB_NAMES.values():
        if os.path.exists(db_name):
            os.remove(db_name)
            print(f"Removed database file: {db_name}")

    # Remove the adsites file
    if os.path.exists(ADSITES_FILENAME):
        os.remove(ADSITES_FILENAME)
        print(f"Removed adsites file: {ADSITES_FILENAME}")

    # Remove the jobs file
    if os.path.exists(JOBS_FILENAME):
        os.remove(JOBS_FILENAME)

    print("All Xsorb run files cleaned up.")
