"""
Created on Tue 28 Feb 2023
@author: Enrico Pedretti

Main functions to generate the adsorption configurations and set up the calculations

"""

from ase.io import read, write
import numpy as np
import os, sys
import glob
#
from ase.constraints import FixCartesian
from slab import Slab
from molecule import Molecule
from io_utils import launch_jobs
from settings import Settings
from dftcode_specific import override_settings, Calculator, OUT_FILE_PATHS
from filenames import *

import ase_custom


#OK (code agnostic)
def adsorption_configurations(settings : Settings, SAVEFIG : bool = False, VERBOSE : bool = False):
    '''
    Generate all adsorption structures considering the combinations of 
    molecular rotations and adsorption sites.

    Returns a list of  ASE atoms with the configurations, 
    and a list of info (csv format) of each configuration.

    Args:
    - settings: Settings object, containing the filenames of molecule and slab, and all the parameters for the rotations of the molecule and the indentifications of adsorption sites
    - SAVEFIG: save an image of the adsorption sites and of the molecular rotations.
    ''' 

    #Slab import from file
    if VERBOSE: print('Loading slab...')
    slab = Slab(slab_filename=settings.slab_filename, 
                layers_threshold=settings.layers_height, 
                surface_sites_height=settings.surface_height, 
                fixed_layers_slab=settings.fixed_layers_slab, 
                fixed_indices_slab=settings.fixed_indices_slab, 
                fix_slab_xyz=settings.fix_slab_xyz,
                sort_atoms_by_z=settings.sort_atoms_by_z,
                translate_slab_from_below_cell_bottom=settings.translate_slab_from_below_cell_bottom)
    if VERBOSE: print('Slab loaded.')


    #Molecule import from file
    if VERBOSE: print('Loading molecule...')
    mol = Molecule(molecule_filename=settings.molecule_filename,
                   atom_index=settings.selected_atom_index,
                   molecule_axis_atoms=settings.molecule_axis_atoms, 
                   axis_vector=settings.axis_vector, 
                   atoms_subset=settings.mol_subset_atoms, 
                   fixed_indices_mol=settings.fixed_indices_mol, 
                   fix_mol_xyz=settings.fix_mol_xyz)
    if VERBOSE: print('Molecule loaded.')


    #Find adsorption sites and labels (site type and x,y coords.)
    adsites, adsites_labels = slab.find_adsorption_sites(
        **settings.sites_find_args,
        selected_sites=settings.selected_sites,
        save_image=SAVEFIG,
        VERBOSE=VERBOSE)


    #Generate all the configs for the various molecular rotations and a list of labels
    all_mol_configs_ase, rotations_labels = mol.generate_molecule_rotations(
        x_rot_angles=settings.x_rot_angles, 
        y_rot_angles=settings.y_rot_angles,
        z_rot_angles=settings.z_rot_angles, 
        vert_angles_list=settings.vertical_angles,
        save_image=SAVEFIG,
        VERBOSE=VERBOSE)


    #Adsorption of molecule on all adsorption sites for all molecule orientations
    if VERBOSE: print('Generating adsorption structures...') 
    all_mol_on_slab_configs_ase = []
    full_labels = []
    for mol_config, rot_label in zip(all_mol_configs_ase, rotations_labels):
        structures, labels = slab.generate_adsorption_structures(molecule=mol_config, 
                                                adsites=adsites,
                                                z_distance_from_site=settings.screening_atom_distance,
                                                min_z_distance_from_surf=settings.screening_min_distance,
                                                adsites_labels=adsites_labels,
                                                rotation_label=rot_label,
                                                mol_before_slab=settings.mol_before_slab) 
        all_mol_on_slab_configs_ase += structures
        full_labels += labels

    if VERBOSE: print('All slab+adsorbate cells generated.')

    return all_mol_on_slab_configs_ase, full_labels

#OK (code agnostic)
def write_labels_csvfile(full_labels : list, labels_filename : str):
    '''
    Write the labels of the configurations to a csv file.

    Args:
    - full_labels: list of labels for the configurations, in csv format
    - labels_filename: name of the csv file where the data will be written
    '''
    with open(labels_filename, 'w') as csvfile:
        csvfile.write('Label' + ',' + 'xrot' + ',' + 'yrot' + ',' + 'zrot' + ',' + 'site' + ',' + 'x' + ',' + 'y' + ',' + 'z' +'\n')
        for i, label in full_labels:
            csvfile.write(f'{i},{label}\n')

#OK (code agnostic)  
def write_inputs(settings : Settings, 
                 all_mol_on_slab_configs_ase : list, 
                 calc_type : str,
                 OVERRIDE_SETTINGS : bool = True, 
                 INTERACTIVE : bool = False):
    '''
    Writes the input files for all the adsorption configurations.

    Returns a list of the indices associated to the files.

    Args:
    - settings: Settings object, containing the dft code parameters
    - all_mol_on_slab_configs_ase: list of ASE atoms for all the adsorption configurations
    - calc_type: 'SCREENING' or 'RELAX'
    - OVERRIDE_SETTINGS: override some specifc settings (e.g. conv tresholds)
    - INTERACTIVE: interactive mode: ask before overwriting files that are already present
    '''
    
    if INTERACTIVE: print('Writing input files...')

    if OVERRIDE_SETTINGS: override_settings(settings, calc_type) 


    ANSWER_ALL = False
    answer = 'yes'
    
    for i, atoms in enumerate(all_mol_on_slab_configs_ase):
        
        file_label = f'{calc_type.lower()}_{i}'

        corresponding_outfile = OUT_FILE_PATHS[calc_type][settings.program].format(i)
        
        #if in interactive mode, ask before overwriting
        if(INTERACTIVE and os.path.isfile(corresponding_outfile)): #convoluted, but it works
            print(f'{corresponding_outfile} already present, possibly from a running calculation. '+ \
                  ('It will be {0}, as requested.'.format('skipped' if 'n' in answer  else 're-calculated') \
                  if ANSWER_ALL else 'You can decide to re-calculate it or skip it.'))
            while True and not ANSWER_ALL:
                answer = input('Re-calculate? ("y" = yes to this one, "yall" = yes to all, "n" = no to this one, "nall" = no to all): ')
                if answer == 'yes' or answer == 'y' or answer == 'yall' or answer == 'no' or answer == 'n' or answer == 'nall': 
                    if answer == 'yall' or answer == 'nall': ANSWER_ALL = True
                    break
                else: print('Value not recognized. Try again.')
            if answer == 'no' or answer == 'n' or answer == 'nall': continue #skip if user does not want to overwrite
        
        calc = Calculator(settings, file_label, atoms) 
        calc.write_input(atoms)

    if INTERACTIVE: print('All input files written.')

    return np.arange(len(all_mol_on_slab_configs_ase))

#OK (code agnostic)  
def generate(SAVEFIG=False):
    '''
    Generates adsorption configurations and writes the inputs, using the settings for the final relaxations.
    Useful to use the program just as a generator of configurations

    Args:
    - SAVEFIG: save an image of the adsorption sites and of the molecular rotations when generating the configurations
    '''
    
    settings=Settings()

    all_mol_on_slab_configs_ase, full_labels = adsorption_configurations(settings, SAVEFIG)

    write_labels_csvfile(full_labels, labels_filename=labels_filename)

    write_inputs(settings, 
                 all_mol_on_slab_configs_ase,
                 calc_type='SCREENING',
                 OVERRIDE_SETTINGS=False)

#OK (code agnostic)   
def launch_screening(SAVEFIG : bool = False):
    '''
    Generates adsorption configurations, writes inputs and launches calculations for the preliminary screening.

    Args:
    - SAVEFIG: save an image of the adsorption sites and of the molecular rotations when generating the configurations
    '''
    
    settings=Settings()

    all_mol_on_slab_configs_ase, full_labels = adsorption_configurations(settings, SAVEFIG, VERBOSE=True)

    write_labels_csvfile(full_labels, labels_filename=labels_filename)

    indices_list = write_inputs(settings, 
                                 all_mol_on_slab_configs_ase,
                                 calc_type='SCREENING', 
                                 VERBOSE=True)

    launch_jobs(program=settings.program,
                calc_type='SCREENING',
                jobscript=settings.jobscript,
                sbatch_command=settings.sbatch_command,
                indices_list=indices_list)    

#TODO: check later that if not relax_completed. Before we put None in the energy
def final_relax(n_configs: int = None, threshold : float = None, exclude : list= None, indices : list = None, REGENERATE=False, BY_SITE = False):
    
    #check: if all calculations finished, simply skip
    #energies = get_energies(pwo_prefix='relax', VERBOSE=False)
    #if energies and None not in energies:
    #    print('All screening calculations completed. Nothing will be done.')
    #    return    
    
    
    if n_configs is None and threshold is None: n_configs = N_relax_default if not BY_SITE else 1
    
    #Check for duplicates in exclude or in indices
    if indices is not None:
        if (len(set(indices)) < len(indices)):
            print('The indices list contain duplicate elements. Quitting.')
            sys.exit(1)
    if exclude is not None:
        if (len(set(exclude)) < len(exclude)):
            print('The exclude list contain duplicate elements. Quitting.')
            sys.exit(1)        


    settings=Settings()

    files_list = [ file.replace('pwi', 'pwo') if os.path.isfile(file.replace('pwi', 'pwo')) else None for file in natsorted(glob.glob( "screening_*.pwi" ))]
    #this will contain None if the pwi has not the corresponding pwo

    if indices is not None: calcs = indices
    else:
        print('Collecting energies from screening...')

        
        energies = [ get_energy_from_pwo(file, REQUIRE_RELAX_COMPLETED=True) if file is not None else None for file in files_list]
        #this will contain None if there are some pwi for which the corresponding pwo does not exist or has not a final relaxed energy

        
        if None in energies:
            print('Not all screening calculations have provided a final energy. Those for which no pwo exists will be excluded, while for those not completed the last coordinates will be used.')
            while True:
                answer = input('Continue anyway with the ones available? ("y" = yes, "n" = no (quit)): ')
                if answer == 'yes' or answer == 'y' or answer == 'no' or answer == 'n': 
                    break
                else: print('Value not recognized. Try again.')
            if 'n' in answer:
                sys.exit(1)
        print('screening energies collected.')
        #if we decided to continue anyway, set the None cases to an extremely high value so that they do not get picked as lowest energy
        #this is just a trick, put in a more elegant way in the future
        energies = [ en if en is not None else 1e50 for en in energies]

        if exclude is None: exclude = []
        else: print('Configurations {0} will be excluded, as requested'.format(exclude))


        calcs = []

        if(BY_SITE):
            site_labels = []
            #TODO: if not present, generate site_labels.csv
            with open('site_labels.csv', 'r') as f:
                file = f.readlines()

                for line in file:
                    if 'site' in line: continue
                    site_labels.append(int(line.split(',')[4].split(' ')[0]))        
            sites = set(site_labels)
            site_labels = np.array(site_labels) 
            samesite_indices = [np.where(site_labels == site) for site in sites] #config labels associated to same site

            for indices_of_site in samesite_indices:
                energies_site = [energies[j] for j in indices_of_site[0] ]
                config_labels_site = [j for j in indices_of_site[0] ]
               
                sorted_indices = np.argsort(energies_site, kind='stable').tolist()
                if n_configs is not None:                  
                    if(n_configs > len(sorted_indices)):
                        print('Error! The number of screening configurations is lower than the desired number of configurations to fully relax. Try with a lower value of --n (the default is 5)')
                        sys.exit(1)
                    for n, i in enumerate(sorted_indices):
                        if config_labels_site[i] in exclude: continue
                        if n < n_configs:
                            calcs.append(config_labels_site[i])
                elif threshold is not None:
                    e_min = min(energies_site)
                    calcs += [config_labels_site[i] for i in sorted_indices if energies_site[i] - e_min <= threshold and config_labels_site[i] not in exclude]


        else:
            sorted_indices = np.argsort(energies, kind='stable').tolist()
            if n_configs is not None:           
                if(n_configs > len(sorted_indices)):
                    print('Error! The number of screening configurations is lower than the desired number of configurations to fully relax. Try with a lower value of --n (the default is 5)')
                    sys.exit(1)
                for i in sorted_indices:
                    if i in exclude: continue
                    if len(calcs)<n_configs:
                        calcs.append(i)
            elif threshold is not None:
                e_min = min(energies)
                calcs += [i for i in sorted_indices if energies[i] - e_min <= threshold and i not in exclude]


    settings.espresso_settings_dict['CONTROL'].update({'calculation' : 'relax'})
    settings.espresso_settings_dict['CONTROL'].update({'restart_mode' : 'from_scratch'})
    settings.espresso_settings_dict['IONS'].update({'ion_dynamics': settings.ion_dynamics})

    #Slab import from file
    slab = Slab(settings.slab_filename, layers_threshold=settings.layers_height, surface_sites_height=settings.surface_height, 
                fixed_layers_slab=settings.fixed_layers_slab, 
                fixed_indices_slab=settings.fixed_indices_slab, 
                fix_slab_xyz=settings.fix_slab_xyz,
                sort_atoms_by_z=True)
    #Molecule import from file
    mol = Molecule(settings.molecule_filename, settings.molecule_axis_atoms, settings.axis_vector, settings.mol_subset_atoms, 
                   settings.fixed_indices_mol, 
                   settings.fix_mol_xyz)

    if(REGENERATE):
        #Re-generation of structures, placing the molecule very close to the surface to avoid
        #trapping in local physisorption minima
        adsites, adsites_labels = slab.find_adsorption_sites(
            *settings.sites_find_args.values(), 
            save_image=False,
            selected_sites=settings.selected_sites)

        all_mol_configs_ase, configs_labels = mol.generate_molecule_rotations(
            atom_index=settings.selected_atom_index,  
            y_rot_angles=settings.y_rot_angles,
            x_rot_angles=settings.x_rot_angles, 
            z_rot_angles=settings.z_rot_angles, 
            vert_angles_list=settings.vertical_angles,
            distance_from_surf=settings.relax_atom_distance, 
            min_distance=settings.relax_min_distance, 
            save_image=False
            )

        #Adsorption of molecule on all adsorption sites for all molecule orientations
        print('Generating adsorption structures...') 
        all_mol_on_slab_configs_ase = [slab.generate_adsorption_structures(molecule=mol_config, adsites=adsites) for mol_config in all_mol_configs_ase]
        all_mol_on_slab_configs_ase = sum(all_mol_on_slab_configs_ase, []) #flatten the 2D array
        if(False): all_mol_on_slab_configs_ase = adsorb_both_surfaces(all_mol_on_slab_configs_ase) #Replicate molecule on the other side of the slab. NOTE: Currently not working!
        full_labels = [mol_config[0]+site_label+mol_config[1] for mol_config in configs_labels for site_label in adsites_labels]
        print('All slab+adsorbate cells generated.')
    else:
        all_mol_on_slab_configs_ase = [None] *  len(files_list)
        for i, file in enumerate(files_list):
            if file is None: continue
            try:
                r = read(filename=file, results_required=False)
            except AssertionError as e:
                print("Error while reading {0} due to ASE problem in reading calculations with scf NOT terminated. Skipping.".format(file))
                continue
            #set constraints from mol and slab if when we are reading the whole pre-relaxed structure
            c = [constraint for constraint in slab.slab_ase.constraints] + \
                [FixCartesian(constraint.a + slab.natoms, ~constraint.mask) for constraint in mol.mol_ase.constraints]
            #for some obscure reason FixCartesian wants the mask in the form 1=Fix, 0=Free, but then it negates in internally, so we first need to negate the mask to construct a new object
            r.set_constraint(c)
            all_mol_on_slab_configs_ase[i] = r

 

    pwi_names = []

    ANSWER_ALL = False
    answer = 'yes'
    
    for i in calcs: #calcs are the integer indexes of the selected calculations to run

        if all_mol_on_slab_configs_ase[i] == None: continue

        file_prefix = 'relax_'+str(i)
        filename = f'{file_prefix}.pwi'

        if(os.path.isfile(filename.replace('pwi', 'pwo'))): 
            print(filename.replace('pwi', 'pwo')+' already present, possibly from a running calculation. '+('It will be {0}, as requested.'.format('skipped' if 'n' in answer else 're-calculated') \
                  if ANSWER_ALL else 'You can decide to re-calculate it or skip it.'))
            while True and not ANSWER_ALL:
                answer = input('Re-calculate? ("y" = yes to this one, "yall" = yes to all, "n" = no to this one, "nall" = no to all): ')
                if answer == 'yes' or answer == 'y' or answer == 'yall' or answer == 'no' or answer == 'n' or answer == 'nall': 
                    if answer == 'yall' or answer == 'nall': ANSWER_ALL = True
                    break
                else: print('Value not recognized. Try again.')
            if answer == 'no' or answer == 'n' or answer == 'nall': continue #skip if user does not want to overwrite


        pwi_names.append(filename)
        calc = Espresso(pseudopotentials=settings.pseudopotentials, 
                    input_data=settings.espresso_settings_dict,
                    label=file_prefix,
                    kpts= settings.kpoints[1] if 'gamma' not in settings.kpoints[0] else None, koffset=settings.kpoints[2] if 'gamma' not in settings.kpoints[0] else None)    
        calc.write_input(all_mol_on_slab_configs_ase[i])

    #print(pwi_names)
    launch_jobs(jobscript=settings.jobscript, pwi_list=pwi_names, outdirs=relax_outdir, jobname_title='rel')


def saveas(which : str, i_or_f : str, saveas_format : str):

    if i_or_f == 'i':
        pw = 'pwi'
    elif i_or_f == 'f':
        pw = 'pwo'
    else: 
        raise RuntimeError("Wrong arguments: passed '{0} {1}', expected 'screening i/f' or 'relax i/f'".format(which, i_or_f))
    if which != 'screening' and which != 'relax':
        raise RuntimeError("Wrong argument: passed '{0}', expected 'screening' or 'relax'".format(which))

    print('Reading files...')
    pw_list=glob.glob(which+'_*.'+pw)
    configs = [(read(file) if pw == 'pwi' else read(file, results_required=False)) for file in pw_list]
    print('All files read.')

    print("Saving {0} files to {1} format...".format(which, saveas_format))
    folder = saveas_format+'/'+which+'/'
    if not os.path.exists(saveas_format):
        os.mkdir(saveas_format)
    if not os.path.exists(folder):
        os.mkdir(folder)    

    for i, config in enumerate(configs):
            if(saveas_format == 'xyz'):
                ase_custom.write_xyz_custom(folder+pw_list[i].replace(pw, saveas_format), config)
            else:
                write(folder+pw_list[i].replace(pw, saveas_format), config)

    print("Files saved to {0}".format(saveas_format+'/'+which) )
    

