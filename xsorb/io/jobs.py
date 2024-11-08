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
import subprocess
import sys

import xsorb.io.database
from xsorb.io.settings import Settings
from xsorb.dft_codes.definitions import SBATCH_POSTFIX
from xsorb.dft_codes.calculator import edit_files_for_restart
if TYPE_CHECKING:
    from xsorb.io.inputs import WrittenSystem


TEST = False


def launch_jobs(*,program : str,
                calc_type : str,
                jobscript : str,
                sbatch_command : str,
                systems : list[WrittenSystem],
                jobname_prefix : str = ''):
    '''
    Launch the calculations.
    Writes the job ids in the database.

    Args:
    - program: 'espresso', 'vasp' or 'ml'
    - calc_type: 'screening'/'relax'/'ml_opt' or 'isolated'
    - jobscript: path of the jobscript file
    - sbatch_command: command to submit the jobscript (in Slurm it is sbatch)
    - systems: list of WrittenSystem objects containing calc_id and paths
    - jobname_prefix: prefix for the job name

    '''
    main_dir = os.getcwd()

    submitted_jobs = []
    for system in systems:

        j_dir = Path(system.in_file_path).parent
        shutil.copyfile(jobscript, f'{j_dir}/jobscript.sh')

        os.chdir(j_dir)   ####################

        #change job title (only for slumr jobscripts)
        with open('jobscript.sh', 'r',encoding=sys.getfilesystemencoding()) as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if "job-name" in line:
                    prefix = jobname_prefix
                    if jobname_prefix != '': prefix += '_' #pylint: disable=multiple-statements
                    if calc_type != 'isolated':
                        suffix = f'{calc_type[0]}{system.calc_id}'
                    else:
                        suffix = system.calc_id
                    lines[i] = f"{line.split('=')[0]}={prefix}{suffix}\n"
                    break
        with open('jobscript.sh', 'w',encoding=sys.getfilesystemencoding()) as f:
            f.writelines(lines)

        postfix = SBATCH_POSTFIX[program].format(
            in_file=Path(system.in_file_path).name,
            out_file=Path(system.out_file_path).name,
            log_file=Path(system.log_file_path).name,
            main_dir=main_dir)
        launch_string = f"{sbatch_command} jobscript.sh {postfix}"
        if TEST: print(launch_string) #pylint: disable=multiple-statements
        else:
            outstring = subprocess.getoutput(launch_string) #launches the jobscript from j_dir
            print(outstring)
            submitted_jobs.append(int(outstring.split()[-1]))
        os.chdir(main_dir) ####################

    if calc_type not in ('isolated'): #no database for slab/molecule
        xsorb.io.database.Database.add_job_ids(calc_type, [system.calc_id for system in systems], submitted_jobs)
    else:
        with open(".submitted_jobs.txt", "a",encoding=sys.getfilesystemencoding()) as f:
            f.writelines([f'{job}\n' for job in submitted_jobs])


def restart_jobs(calc_type : str):
    '''
    Restart the uncompleted dft calculations.
    Associated to the command 'xsorb restart screening/relax' in the CLI.
    Beware:n o restart for ML!

    Args:
    - calc_type: 'screening' or 'relax'.
    '''

    settings = Settings()

    rows = xsorb.io.database.Database.get_calculations(calc_type,
                                     selection='status="incomplete",job_status="terminated"')
    indices_to_restart = [row.calc_id for row in rows]
    in_files = [row.in_file for row in rows]
    out_files = [row.out_file for row in rows]
    log_files = [row.log_file for row in rows]

    #edit input files
    edit_files_for_restart(settings.program, in_files)

    #launch the calculations
    main_dir = os.getcwd()
    submitted_jobs = []
    for in_file, out_file, log_file in zip(in_files, out_files, log_files):

        j_dir = Path(in_file).parent
        os.chdir(j_dir)

        postfix = SBATCH_POSTFIX[settings.program].format(in_file=in_file,
                                                 out_file=out_file,
                                                 log_file=log_file,
                                                 main_dir=main_dir)
        launch_string = f"{settings.input.submit_command} jobscript.sh {postfix}"

        if TEST: print(launch_string) #pylint: disable=multiple-statements
        else:
            outstring = subprocess.getoutput(launch_string) #launches the jobscript from j_dir
            print(outstring)
            submitted_jobs.append(int(outstring.split()[-1]))
        os.chdir(main_dir) ####################

    xsorb.io.database.Database.add_job_ids(calc_type, indices_to_restart, submitted_jobs)


#TODO: generalize to other schedulers, replace with regex
def get_running_jobs():
    '''
    Get the running jobs by interrogating the scheduler
    '''
    running_jobs = subprocess.getoutput("squeue --me").split("\n")[1:]
    running_job_ids = [int(job.split()[0]) for job in running_jobs]

    return running_job_ids


def _cancel_jobs(ids : list[int]):
    '''
    Cancel the jobs with the given ids
    '''
    for job_id in ids:
        os.system(f"scancel {job_id}")


def scancel():
    '''
    Cancel all the running jobs For the current Xsorb session.
    Associated to the command 'xsorb scancel' in the CLI.
    '''

    #add jobs from the database(s)
    submitted_job_ids = xsorb.io.database.Database.get_all_job_ids()

    #also add jobs from .submitted_jobs.txt
    if Path(".submitted_jobs.txt").exists():
        with open(".submitted_jobs.txt", "r",encoding=sys.getfilesystemencoding()) as f:
            submitted_jobs = f.readlines()
            submitted_job_ids.extend([job.strip() for job in submitted_jobs])

    running_jobs = get_running_jobs()
    job_ids_to_cancel = [job for job in running_jobs if job in submitted_job_ids]

    if len(job_ids_to_cancel) == 0:
        print("No jobs to cancel.")
        return

    print(f"Cancelling jobs {job_ids_to_cancel}.")
    _cancel_jobs(job_ids_to_cancel)
    print("All jobs cancelled.")
