'''
Moudule to clean all files for an Xsorb run.
'''

import os
import shutil
import logging

from xsorb.io.filenames import (
    ALL_DB_NAMES, JOBS_FILENAME,
    ALL_OUTDIRS, ALL_FILES_AND_DIRS)
from xsorb.io.utils import yes_no_question

def cleanup_xsorb_run(calc_only: bool = False):
    '''
    Clean all files for an Xsorb run.
    '''

    # ask the user for confirmation
    if not yes_no_question("Do you really want to remove all files for this Xsorb run?"):
        logging.info("Cleanup aborted.")
        return


    logging.info("Cleaning up calculation files...")
    for directory in ALL_OUTDIRS:
        if os.path.exists(directory):
            shutil.rmtree(directory, ignore_errors=True)
            logging.info(f"Removed directory: {directory}")
    logging.info("Calculation files cleaned up.")
    if calc_only:
        return

    # Remove all database files
    for db_name in ALL_DB_NAMES:
        if os.path.exists(db_name):
            os.remove(db_name)
            logging.info(f"Removed database file: {db_name}")

    # Remove the jobs file
    if os.path.exists(JOBS_FILENAME):
        os.remove(JOBS_FILENAME)

    logging.info("All Xsorb run files cleaned up.")


def fresh_start():
    '''
    Check if any of the files from a previous run exists,
    and ask the user if they want to continue or clean up.
    '''

    if any(os.path.exists(file_or_dir) for file_or_dir in ALL_FILES_AND_DIRS):
        if yes_no_question('''Warning: some files from a previous run exist.
              Do you want to delete all of them and start a fresh run?'''):
            cleanup_xsorb_run()
