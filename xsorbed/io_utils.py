#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue 28 Feb 2023
@author: Enrico Pedretti

Function definitions to read from pwo and launch scripts

"""

import os, sys, shutil
import glob
from natsort import natsorted

from xsorbed.filenames import *

TEST = False #do not actually launch the jobs, simply prints the command



def get_energies(in_filename : str, out_filename : str, E_slab_mol : list, pwo_prefix : str):
    #can be called before all the jobs have finished

    rydbergtoev = 13.605703976

    #Begin script
    files = natsorted(glob.glob( pwo_prefix + "_*.pwo" ))

    with open(in_filename, 'r') as f:
        data = f.readlines()


    #add E_ads column header
    line = data[0].split(',')
    line[-1] = line[-1].split('\n')[0]
    if "relax" in pwo_prefix: #relax case
        if E_slab_mol:
            line.append('Eads_rel(eV)\n')
        else:
            line.append('Etot_rel(eV)\n')
    else: #scf case
        if E_slab_mol:
            line.append('Eads_scf(eV)\n')
        else:
            line.append('Etot_scf(eV)\n')
    data[0] = ','.join(line)


    #get energies from pwo(s)
    energies = len(files)*[None]  #those not completed will be left as None, raising an error if used for relax

    for i, file in enumerate(files):
        
        with open(file, 'r') as f:
            pwo = f.readlines()

            scf_terminated = False
            relax_terminated = False
            for line in pwo: #make sure to get the last one (useful in relaxations)
                if '!' in line: 
                    toten = line.split()[4]
                if 'convergence has been achieved' in line:
                    scf_terminated = True
                if 'Final energy' in line:
                    relax_terminated = True


            if scf_terminated: #add energy to line in csv

                config_label = int( (file.split('.pwo')[0]).split('_')[-1] )                

                toten = float(toten)
                if E_slab_mol:
                    toten -= (E_slab_mol[0]+E_slab_mol[1])
                
                toten *= rydbergtoev
                energies[i] = toten


                line = data[config_label+1].split(',')
                line[-1] = line[-1].split('\n')[0]
                line.append('{:.3f}'.format(toten))
                if "relax" in pwo_prefix and not relax_terminated:
                    print(file.split('/')[-1] + ' relaxation has not reached final configuration. The energy will be marked with a *')
                    line[-1]+='*'
                
                data[config_label+1] = ','.join(line)
                data[config_label+1] += '\n'
                
            else: 
                print(file.split('/')[-1] + ' job has not reached scf convergence. It will be skipped.')

    with open(out_filename, 'w') as f:
        f.writelines( data )

    return energies


def get_z(pwo_filename : str, atom_index : int):

    with open(pwo_filename, 'r') as f:

        pwo = f.readlines()

        relax_terminated = False
        first_index = 0
        for i, line in enumerate(pwo):
            if 'Final energy' in line:
                relax_terminated = True
                first_index  = i+3
                #do not break, so that if more than one relax is present we take the last one

        if not relax_terminated: raise RuntimeError(pwo_filename + ' relax not terminated. Quitting.')
        else: 
            z = pwo[first_index+atom_index].split()[3]           
            return float(z)
            

def launch_jobs(jobscript : str, pwi_list : list, outdirs : str, jobname_prefix : str, pwi_prefix : str, pwo_prefix : str):
    #NOTE: for this script, the execution of pw.x in the jobscript need to be called with this option:  
    #      -input $1 >> $2, so that it reads the input from file $1 and writes output in file $2    

    main_dir = os.getcwd()
    if not os.path.exists(outdirs): os.mkdir(outdirs)

    for input_file in pwi_list:

        output_file = input_file.replace("pwi", "pwo")
        output_file = output_file.replace(pwi_prefix, pwo_prefix)
        label = input_file.split('.pwi')[0].split('_')[-1]

        if(os.path.isfile(input_file)): # unnecessary, this check is also done before calling the function

            if(os.path.isfile(output_file)): 
                print(output_file+' already present, possibly from a running calculation. It will be skipped.')
                continue

            j_dir = outdirs+'/'+str(label)
            os.mkdir(j_dir)

            shutil.copyfile(jobscript, j_dir+'/'+jobscript_filename)

            os.chdir(j_dir) #####################
            with open(jobscript_filename, 'r') as f:
                lines = f.readlines()

                for i, line in enumerate(lines):
                    if "job-name" in line:
                        lines[i] = line.split('=')[0] + '="' + jobname_prefix + '_' + label +'"\n'
                        break
        
            with open(jobscript_filename, 'w') as f:
                f.writelines( lines )


            if(TEST): print("sbatch " + jobscript_filename + ' ' +main_dir+'/'+input_file + ' '+ main_dir+'/'+output_file)
            else: os.system("sbatch " + jobscript_filename + ' ' +main_dir+'/'+input_file + ' '+ main_dir+'/'+output_file)  #launchs the jobscript in j_dir from j_dir
            os.chdir(main_dir) ####################


def _is_completed(pwo : str, which : str):
    if(which == 'scf'):
        searchfor = 'End of self-consistent calculation'
    elif(which == 'relax' or which == 'prerelax'):
        searchfor = 'Final energy'
    
    with open(pwo, 'r') as f:
        file = f.readlines()

    completed = False
    for line in file:
        if searchfor in line:
            completed = True
    
    return completed


def restart_jobs(which : str, pwi_prefix : str, pwo_prefix : str):
    if(which == 'scf'):
        outdirs = 'scf_outdirs'
        pwo_prefix_full = pwo_prefix + 'scf'
    elif(which == 'relax'):
        outdirs = 'relax_outdirs'
        pwo_prefix_full = pwo_prefix + 'relax'
    else:
        raise ValueError("Not clear which calculation should be restarted.")
    

    main_dir = os.getcwd()
    pwos = natsorted(glob.glob( pwo_prefix_full + "_*.pwo" ))
    #restart only non-completed calculations
    pwos = [pwo for pwo in pwos if not _is_completed(pwo, which=which)]
    pwis = [pwo.replace(pwo_prefix, pwi_prefix).replace('.pwo', '.pwi') for pwo in pwos]
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
        if(TEST): print("sbatch " + jobscript_filename + ' ' +main_dir+'/'+pwis[i] + ' '+ main_dir+'/'+pwos[i])
        else: os.system("sbatch " + jobscript_filename + ' ' +main_dir+'/'+pwis[i] + ' '+ main_dir+'/'+pwos[i])  #launchs the jobscript in j_dir from j_dir
        os.chdir(main_dir) ####################