#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue 28 Feb 2023
@author: Enrico Pedretti

Function definitions to read from pwo and launch scripts

"""

import os, shutil
import glob
from natsort import natsorted

TEST = False #do not actually launch the jobs, simply prints the command


def get_energies(labels_filename : str, energies_filename : str, pwo_prefix : str = 'output'):
    #can be called before all the jobs have finished

    rydbergtoev = 13.605703976

    #Begin script
    files = natsorted(glob.glob( pwo_prefix + "_*.pwo" ))

    with open(labels_filename, 'r') as f:
        data = f.readlines()

    #add E_ads column header
    line = data[0].split(',')
    line[-1] = line[-1].split('\n')[0]
    if "final" in pwo_prefix:
        line.append('E_rel(eV)\n')
    else:
        line.append('E_scf(eV)\n')
    data[0] = ','.join(line)

    energies = len(files)*[None]

    for i, file in enumerate(files):
        
        with open(file, 'r') as f:
            pwo = f.readlines()

            job_finished = False
            scf_terminated = False
            for line in pwo: #make sure to get the last one (useful in relaxations)
                if '!' in line: 
                    toten = line.split()[4]
                if 'JOB DONE.' in line:
                    job_finished = True
                if 'End of self-consistent calculation' in line:
                    scf_terminated = True


            if job_finished and scf_terminated: #add energy to line in csv

                config_label = int( (file.split('.pwo')[0]).split('_')[-1] )
                
                toten = float(toten)*rydbergtoev
                energies[i] = toten

                line = data[config_label+1].split(',')
                line[-1] = line[-1].split('\n')[0]
                line.append('{:.3f}'.format(toten))
                
                data[config_label+1] = ','.join(line)
                data[config_label+1] += '\n'
            else: 
                print(file.split('/')[-1] + ' job is not (yet) completed. It will be skipped.')


    with open(energies_filename, 'w') as f:
        f.writelines( data )

    return energies


def get_z(pwo_filename : str, atom_index : int):

    with open(pwo_filename, 'r') as f:

        pwo = f.readlines()

        job_finished = False
        first_index = 0
        for i, line in enumerate(pwo): #make sure to get the last one (useful in relaxations)
            if 'Begin final coordinates' in line:
                job_finished = True
                first_index  = i+3
                #do not break, so that if more than one run is present we take the last one

        if not job_finished: raise RuntimeError(pwo_filename + ' job did not terminate correctly. Quitting.')
        else: 
            z = pwo[first_index+atom_index].split()[3]           
            return float(z)
            

def launch_jobs(jobscript : str, pwi_list : list[str], outdirs : str, jobname_prefix=''):
    #NOTE: for this script, the execution of pw.x in the jobscript need to be called with this option:  
    #      -input $1 >> $2, so that it reads the input from file $1 and writes output in file $2    

    main_dir = os.getcwd()
    if not os.path.exists(outdirs): os.mkdir(outdirs)

    for input_file in pwi_list:

        output_file = input_file.replace("pwi", "pwo")
        output_file = output_file.replace("input", "output")
        label = input_file.split('.pwi')[0].split('_')[-1]

        if(os.path.isfile(input_file)): # unnecessary, this check is also done before calling the function

            if(os.path.isfile(output_file)): 
                raise RuntimeError(output_file+' already present, possibly from a running calculation. Quitting.')

            j_dir = outdirs+'/'+str(label)
            os.mkdir(j_dir)

            shutil.copyfile(jobscript, j_dir+'/jobscript')

            os.chdir(j_dir) #####################
            with open("jobscript", 'r') as f:
                lines = f.readlines()

                for i, line in enumerate(lines):
                    if "job-name" in line:
                        lines[i] = line.split('=')[0] + '="' + jobname_prefix + '_' + label +'"\n'
                        break
                

            with open("jobscript", 'w') as f:
                f.writelines( lines )

            if(TEST): print("sbatch jobscript" + ' ' +main_dir+'/'+input_file + ' '+ main_dir+'/'+output_file)
            else: os.system("sbatch jobscript " + ' ' +main_dir+'/'+input_file + ' '+ main_dir+'/'+output_file)  #launchs the jobscript in j_dir from j_dir
            os.chdir(main_dir) ####################


def restart_jobs(which : str = 'scf'):
    if(which == 'scf'):
        outdirs = 'scf_outdirs'
        pwo_prefix = 'output_scf'
    elif(which == 'prerelax'):
        outdirs = 'prerel_outdirs'
        pwo_prefix = 'output_prerelax'
    elif(which == 'finalrelax'):
        outdirs = 'finalrel_outdirs'
        pwo_prefix = 'output_finalrelax'
    else:
        raise ValueError("Not clear which calculation should be restarted.")
    

    main_dir = os.getcwd()
    pwos = natsorted(glob.glob( pwo_prefix + "_*.pwo" ))
    pwis = [pwo.replace('output', 'input').replace('.pwo', '.pwi') for pwo in pwos]
    config_labels = [(file.split('.pwo')[0]).split('_')[-1] for file in pwos]


    #edit pwi(s)
    for pwi in pwis:
        with open(pwi, 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if 'from_scratch' in line:
                    lines[i] = lines[i].replace('from_scratch','restart')
                    break
        with open(pwi, 'w') as f:
            f.writelines( lines )    

    #launch jobs
    for i, label in enumerate(config_labels):        
        os.chdir(outdirs+'/'+label) #####################
        if(TEST): print("sbatch jobscript" + ' ' +main_dir+'/'+pwis[i] + ' '+ main_dir+'/'+pwos[i])
        else: os.system("sbatch jobscript " + ' ' +main_dir+'/'+pwis[i] + ' '+ main_dir+'/'+pwos[i])  #launchs the jobscript in j_dir from j_dir
        os.chdir(main_dir) ####################