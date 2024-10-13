'''
Module for launching the calculations
'''

import os
from pathlib import Path
import shutil
import subprocess

from xsorb import common_definitions
from xsorb.settings import Settings
from xsorb.dft_codes.definitions import SBATCH_POSTFIX,IN_FILE_PATHS,OUT_FILE_PATHS,LOG_FILE_PATHS
from xsorb.dft_codes.calculator import edit_files_for_restart
from xsorb.io.database import Database


TEST = True


calc_dirs = {'screening': common_definitions.screening_outdir,
             'relax': common_definitions.relax_outdir,
             'ml_opt': common_definitions.preopt_outdir}


#TODO: generalize to other schedulers, replace with regex
def get_running_jobs():
    '''
    Get the running jobs by interrogating the scheduler
    '''
    running_jobs = subprocess.getoutput("squeue --me").split("\n")[1:]
    running_job_ids = [job.split()[0] for job in running_jobs]

    return running_job_ids

def cancel_jobs(ids : list[int]):
    '''
    Cancel the jobs with the given ids
    '''
    for job_id in ids:
        os.system(f"scancel {job_id}")


def launch_jobs(program : str,
                calc_type : str,
                jobscript : str,
                sbatch_command : str,
                systems : list[dict],
                jobname_prefix : str = ''):
    '''
    Launch the calculations.
    Writes the job ids in the database.

    Args:
    - program: 'espresso', 'vasp' or 'ml'
    - calc_type: 'screening', 'relax' or 'preopt'
    - jobscript: path of the jobscript file
    - sbatch_command: command to submit the jobscript (in Slurm it is sbatch)
    - systems: list of dictionaries containing the calc_id key
    - jobname_prefix: prefix for the job name

    '''
    main_dir = os.getcwd()

    submitted_jobs = []
    for system in systems:

        calc_id = system['calc_id']

        j_dir = f'{calc_dirs[calc_type]}/{calc_id}'

        shutil.copyfile(jobscript, j_dir)

        os.chdir(j_dir)   ####################

        #change job title (only for slumr jobscripts)
        with open(Path(jobscript).name, 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if "job-name" in line:
                    prefix = jobname_prefix
                    if jobname_prefix != '': prefix += '_' #pylint: disable=multiple-statements
                    lines[i] = f"{line.split('=')[0]}={prefix}{calc_type[0]}{calc_id}\n"
                    break
        with open(Path(jobscript).name, 'w') as f:
            f.writelines(lines)

        postfix = SBATCH_POSTFIX[program].format(
            in_file=Path(IN_FILE_PATHS[calc_type][program].format(calc_id)).name,
            out_file=Path(OUT_FILE_PATHS[calc_type][program].format(calc_id)).name,
            log_file=Path(LOG_FILE_PATHS[calc_type][program].format(calc_id)).name,
            main_dir=main_dir)
        launch_string = f"{sbatch_command} {Path(jobscript).name} {postfix}"
        if TEST: print(launch_string) #pylint: disable=multiple-statements
        else:
            outstring = subprocess.getoutput(launch_string) #launches the jobscript from j_dir
            print(outstring)
            submitted_jobs.append(outstring.split()[-1])
        os.chdir(main_dir) ####################

    Database.add_job_ids(calc_type, [system['calc_id'] for system in systems], submitted_jobs)


def launch_independend_jobs(program : str,
                            jobscript : str,
                            sbatch_command : str,
                            systems : list[dict],
                            jobname_prefix : str = ''):
    '''
    Launch custom calculations, outside of screening/relax,
    without touching the database. Writes the job ids in a .submitted_jobs.txt file.

    Args:
    - program: DFT program. Possible values: 'ESPRESSO' or 'VASP'
    - jobscript: path of the jobscript file
    - sbatch_command: command to submit the jobscript (in Slurm it is sbatch)
    - systems: list of dictionaries containing 'calc_dir', 'label'
    '''
    main_dir = os.getcwd()

    submitted_jobs = []
    for system in systems:

        j_dir = system['calc_dir']
        label = system['label']

        shutil.copyfile(jobscript, j_dir)

        os.chdir(j_dir)   ####################

        #change job title (only for slumr jobscripts)
        with open(Path(jobscript).name, 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if "job-name" in line:
                    prefix = jobname_prefix
                    if jobname_prefix != '': prefix += '_' #pylint: disable=multiple-statements
                    lines[i] = f"{line.split('=')[0]}={prefix}{label[:4]}\n"
                    break
        with open(Path(jobscript).name, 'w') as f:
            f.writelines(lines)

        postfix = SBATCH_POSTFIX[program].format(
            in_file=label+'.xyz',
            out_file=label+'.traj',
            log_file=label+'.log',
            main_dir=main_dir)
        launch_string = f"{sbatch_command} {Path(jobscript).name} {postfix}"
        if TEST: print(launch_string) #pylint: disable=multiple-statements
        else:
            outstring = subprocess.getoutput(launch_string) #launches the jobscript from j_dir
            print(outstring)
            submitted_jobs.append(outstring.split()[-1])
        os.chdir(main_dir) ####################

    with open(".submitted_jobs.txt", "a") as f:
        f.writelines([job+'\n' for job in submitted_jobs])


def restart_jobs(calc_type : str):
    '''
    Restart the uncompleted calculations

    Args:
    - calc_type: 'screening' or 'relax'
    '''

    Database.update_calculations(calc_type)
    indices_to_restart = [row.calc_id for row in \
                          Database.get_calculations(
                              calc_type,
                              selection='status="incomplete",job_status="terminated"')]

    settings = Settings(verbose=False)
    program = settings.program
    main_dir = os.getcwd()

    #edit input files
    edit_files_for_restart(settings.program, calc_type, indices_to_restart)

    #launch the calculations
    submitted_jobs = []
    for calc_id in indices_to_restart:
        j_dir = f'{calc_dirs[calc_type]}/{calc_id}'
        os.chdir(j_dir)

        postfix = SBATCH_POSTFIX[program].format(
            in_file=Path(IN_FILE_PATHS[calc_type][program].format(calc_id)).name,
            out_file=Path(OUT_FILE_PATHS[calc_type][program].format(calc_id)).name,
            log_file=Path(LOG_FILE_PATHS[calc_type][program].format(calc_id)).name,
            main_dir=main_dir)
        launch_string = f"{settings.input.submit_command} "\
            f"{Path(settings.input.jobscript_path).name} {postfix}"
        if TEST: print(launch_string) #pylint: disable=multiple-statements
        else:
            outstring = subprocess.getoutput(launch_string) #launches the jobscript from j_dir
            print(outstring)
            submitted_jobs.append(outstring.split()[-1])
        os.chdir(main_dir) ####################

    Database.add_job_ids(calc_type, indices_to_restart, submitted_jobs)


def scancel_jobs():
    '''
    Cancel all the jobs For the current Xsorb session.
    '''

    #add jobs from the database(s)
    submitted_job_ids = Database.get_all_job_ids()

    #also add jobs from .submitted_jobs.txt
    with open(".submitted_jobs.txt", "r") as f:
        submitted_jobs = f.readlines()
        submitted_job_ids.extend([job.strip() for job in submitted_jobs])

    running_jobs = get_running_jobs()
    running_job_ids = [job.split()[0] for job in running_jobs]

    print(f"Cancelling jobs {running_job_ids}.")
    cancel_jobs(running_job_ids)
    print("All jobs cancelled.")
