#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Enrico Pedretti

'''
Module for launching the calculations
'''

from __future__ import annotations
from typing import TYPE_CHECKING
import os
from pathlib import Path
import shutil
import sys
import logging

import xsorb.io.database
from xsorb.settings import Settings
from xsorb.dft_codes.definitions import SBATCH_POSTFIX
from xsorb.dft_codes.calculator import edit_files_for_restart
from xsorb.io.scheduler import JobScheduler
from xsorb.io.filenames import JOBS_FILENAME
if TYPE_CHECKING:
    from xsorb.adsorptiondata.adsorptioncalculation import AdsorptionCalculation


TEST = False


def launch_jobs(*,program : str,
                calc_type : str,
                jobscript : str,
                scheduler_name : str,
                systems : list[AdsorptionCalculation],
                jobname_prefix : str = ''):
    '''
    Launch the calculations.
    Writes the job ids in the database

    Args:
    - program: 'espresso', 'vasp' or 'ml'
    - calc_type: 'screening'/'relax'/'mlopt' or 'isolated'
    - jobscript: path of the jobscript file
    - scheduler_name: name of the scheduler, e.g. 'slurm'
    - systems: list of WrittenSystem objects containing calc_id and paths
    - jobname_prefix: prefix for the job name

    '''

    scheduler = JobScheduler(scheduler_name)

    main_dir = os.getcwd()

    submitted_jobs : list[str] = []
    for system in systems:

        j_dir = Path(system.calc_info.in_file_path).parent
        shutil.copyfile(jobscript, f'{j_dir}/jobscript.sh')

        os.chdir(j_dir)   ####################

        #change job title (only for slumr jobscripts)
        if scheduler.scheduler_name == 'slurm':
            with open('jobscript.sh', 'r',encoding=sys.getfilesystemencoding()) as f:
                lines = f.readlines()
                for i, line in enumerate(lines):
                    if "job-name" in line:
                        prefix = jobname_prefix[:4]
                        if jobname_prefix != '': prefix += '_' #pylint: disable=multiple-statements
                        if calc_type != 'isolated':
                            suffix = f'{calc_type[0]}{system.calc_info.calc_id}'
                        else:
                            suffix = system.calc_info.calc_id
                        lines[i] = f"{line.split('=')[0]}={prefix}{suffix}\n"
                        break
            with open('jobscript.sh', 'w',encoding=sys.getfilesystemencoding()) as f:
                f.writelines(lines)

        postfix = SBATCH_POSTFIX[program].format(
            in_file=Path(system.calc_info.in_file_path).name,
            out_file=Path(system.calc_info.out_file_path).name,
            log_file=Path(system.calc_info.log_file_path).name,
            main_dir=main_dir)

        jobid = scheduler.submit_job(script_path='jobscript.sh', script_args=postfix.split())
        submitted_jobs.append(jobid)

        os.chdir(main_dir) ####################

    if calc_type not in ('isolated'): #no database for slab/molecule
        xsorb.io.database.Database.add_job_ids(calc_type,
                                               [int(system.calc_info.calc_id) for system in systems],
                                               submitted_jobs)
    else:
        with open(JOBS_FILENAME, "a",encoding=sys.getfilesystemencoding()) as f:
            f.writelines([f'{job}\n' for job in submitted_jobs])


def restart_jobs(calc_type : str):
    '''
    Restart the uncompleted dft calculations.
    Associated to the command 'xsorb restart screening/relax' in the CLI.
    Beware:no restart for ML!

    Args:
    - calc_type: 'screening' or 'relax'.
    '''

    settings = Settings(verbose=False)
    scheduler = JobScheduler(settings.input.scheduler)
    active_jobs = scheduler.get_active_job_ids()

    rows = xsorb.io.database.Database.get_calculations(calc_type=calc_type,
                                     selection='status!=completed')
    indices_to_restart = [row.calc_id for row in rows if row.job_id not in active_jobs]
    in_files = [row.in_file_path for row in rows]
    out_files = [row.out_file_path for row in rows]
    log_files = [row.log_file_path for row in rows]

    #edit input files
    edit_files_for_restart(settings.dft.program, in_files)

    #launch the calculations
    main_dir = os.getcwd()
    submitted_jobs : list[str] = []
    for in_file, out_file, log_file in zip(in_files, out_files, log_files):

        j_dir = Path(in_file).parent
        os.chdir(j_dir)

        postfix = SBATCH_POSTFIX[settings.dft.program].format(in_file=Path(in_file).name,
                                                 out_file=Path(out_file).name,
                                                 log_file=Path(log_file).name,
                                                 main_dir=main_dir)

        jobid = scheduler.submit_job(script_path='jobscript.sh', script_args=postfix.split())
        submitted_jobs.append(jobid)

        os.chdir(main_dir) ####################

    xsorb.io.database.Database.add_job_ids(calc_type, indices_to_restart, submitted_jobs)


def cancel_jobs():
    '''
    Cancel all the running jobs For the current Xsorb session.
    Associated to the command 'xsorb cancel' in the CLI.
    '''

    settings = Settings(verbose=False)
    scheduler = JobScheduler(settings.input.scheduler)

    #add jobs from the database(s)
    submitted_job_ids = xsorb.io.database.Database.get_all_job_ids()

    #also add jobs from .submitted_jobs.txt
    if Path(JOBS_FILENAME).exists():
        with open(JOBS_FILENAME, "r",encoding=sys.getfilesystemencoding()) as f:
            submitted_jobs = f.readlines()
            submitted_job_ids.extend([int(job.strip()) for job in submitted_jobs])

    active_jobs = scheduler.get_active_job_ids()
    job_ids_to_cancel = [job for job in active_jobs if job in submitted_job_ids]

    if len(job_ids_to_cancel) == 0:
        logging.info("No jobs to cancel.")
        return

    logging.info(f"Cancelling jobs {job_ids_to_cancel}.")

    for job_id in job_ids_to_cancel:
        scheduler.cancel_job(job_id)

    logging.info("All jobs cancelled.")
