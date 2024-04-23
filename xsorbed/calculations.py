"""
Created on Tue 28 Feb 2023
@author: Enrico Pedretti

Main functions to generate the adsorption configurations and set up the calculations

"""

from ase.io import read
import numpy as np
import os, sys
from operator import itemgetter
#
from ase.constraints import FixCartesian
from xsorbed.slab import Slab
from xsorbed.molecule import Molecule
from xsorbed.io_utils import launch_jobs, get_calculations_results, _get_configurations_numbers
from xsorbed.settings import Settings
from xsorbed.dftcode_specific import override_settings, Calculator, OUT_FILE_PATHS
from xsorbed.common_definitions import *

from xsorbed import ase_custom



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
                translate_slab_from_below_cell_bottom=settings.translate_slab)
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
        individual_rotations=settings.individual_rotations,
        save_image=SAVEFIG,
        VERBOSE=VERBOSE)

    
    # Move possible additional rotations to the end, to keep the already present calculations untouched
    # BEWARE: removing a rotation breaks the already present indices, possible future TODO: use a dictionary to keep track of the indices
    if os.path.isfile(labels_filename):    
        previous_labels = np.genfromtxt(labels_filename, delimiter=',', names=True)
        previous_rotations = np.unique([x.tolist()[1:4] for x in previous_labels], axis=0)
        new_rotations = [list(map(float, l.split(',')[:-1])) for l in rotations_labels]
        

        #move the new rotations to the end of the list
        for newrot, rot_label_csv, mol_config in zip(new_rotations, rotations_labels.copy(), all_mol_configs_ase.copy()):
        
            if not np.any([np.allclose(newrot, oldrot, atol=1e-3) for oldrot in previous_rotations]):
                
                print(f'New rotation {newrot} found. It will be added to the configurations.')
                rotations_labels.remove(rot_label_csv)
                rotations_labels.append(rot_label_csv) 
                all_mol_configs_ase.remove(mol_config)
                all_mol_configs_ase.append(mol_config)             


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


def write_labels_csvfile(full_labels : list, labels_filename : str):
    '''
    Write the labels of the configurations to a csv file.

    Args:
    - full_labels: list of labels for the configurations, in csv format
    - labels_filename: name of the csv file where the data will be written
    '''
    with open(labels_filename, 'w') as csvfile:
        csvfile.write('Label' + ',' + 'xrot' + ',' + 'yrot' + ',' + 'zrot' + ',' + 'site' + ',' + 'x' + ',' + 'y' + ',' + 'z' +'\n')
        for i, label in enumerate(full_labels):
            csvfile.write(f'{i},{label}\n')


def regenerate_missing_sitelabels():
    '''
    Re-generate the site_labels.csv file if it was accidentally deleted, using the info from settings.in
    '''
    
    settings=Settings()
    
    _, full_labels = adsorption_configurations(settings, SAVEFIG=False, VERBOSE=True)

    write_labels_csvfile(full_labels, labels_filename=labels_filename)

 
def write_inputs(settings : Settings, 
                 all_mol_on_slab_configs_ase : list, 
                 calc_indices : list,
                 calc_type : str,
                 OVERRIDE_SETTINGS : bool = True, 
                 INTERACTIVE : bool = False):
    '''
    Writes the input files for all the adsorption configurations.


    Args:
    - settings: Settings object, containing the dft code parameters
    - all_mol_on_slab_configs_ase: list of ASE atoms for all the adsorption configurations
    - calc_type: 'SCREENING' or 'RELAX'
    - OVERRIDE_SETTINGS: override some specifc settings (e.g. conv tresholds)
    - INTERACTIVE: interactive mode: ask before overwriting files that are already present

    Returns:
    - written_indices: list of the indices of the configurations for which the input files have been written
    '''
    
    if INTERACTIVE: print('Writing input files...')

    if OVERRIDE_SETTINGS: override_settings(settings, calc_type) 


    ANSWER_ALL = False
    answer = 'yes'
    
    written_indices = []

    for i, atoms in zip(calc_indices, all_mol_on_slab_configs_ase):
        
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
        
        
        j_dir = f'{screening_outdir if calc_type == "SCREENING" else relax_outdir}/{i}'
        calc = Calculator(settings, file_label, atoms, j_dir) 
        calc.write_input(atoms)
        written_indices.append(i)

    if INTERACTIVE: print('All input files written.') 

    return written_indices


def obtain_fullrelax_indices(settings : Settings,
                             n_configs: int = None, 
                             threshold : float = None, 
                             exclude : list= None, 
                             BY_SITE = False):
    '''
    Returns a list with the indices of the configurations for the final relaxation,
    chosen according to the given parameters

    Args:
    - settings: Settings object, containing the program type
    - n_configs: nubmer of configurations to be relaxed, starting from the one with lowest energy
    - threshold: energy threshold (in eV) from the NOT EXCLUDED lowest energy configuration. The configuration with E - Emin < threshold will be selected
    - exclude: indices of the configurations to be excluded
    - BY_SITE: do the configuration identification separately for each site. One or more configuration for each site will be produced
    ''' 
    
        
    print('Collecting energies from screening...')
    screening_results = get_calculations_results(program=settings.program, calc_type='SCREENING')
    print('Screening energies collected.')
    
    
    #Ask to quit if some screening calculations are missing
    if False in screening_results['relax_completed'].values() \
        or len(screening_results['relax_completed'].values()) < len(_get_configurations_numbers()):
        
        print('Not all screening calculations have provided a final energy. \
                Those for which no output exists will be excluded, while for those not completed the last coordinates will be used.')
        while True:
            answer = input('Continue anyway with those available? ("y" = yes, "n" = no (quit)): ')
            if answer == 'yes' or answer == 'y' or answer == 'no' or answer == 'n': 
                break
            else: print('Value not recognized. Try again.')
        if 'n' in answer:
            print('Quitting.')
            sys.exit(0)


    if exclude is None: exclude = []
    else: 
        exclude = set(exclude)
        print(f'Configurations {exclude} will be excluded, as requested.')


    #find the indices of the final relaxations
    if(BY_SITE):
        
        site_labels = []
        with open(labels_filename, 'r') as f:
            file = f.readlines()[1:]
            for line in file:
                site_labels.append(int(line.split(',')[4].split(' ')[0]))
        site_labels = np.array(site_labels)       

        #list of lists, outer: site label, inner: calc indices at that site
        site_grouped_index_list = [np.where(site_labels == site)[0] for site in set(site_labels)] 

        
        calculations_indices = []
        for indiceslist_at_site in site_grouped_index_list:                
            #list of {index: energy} for the indices at the site
            indices_and_energies_at_site = [(index, energy) for index, energy in screening_results['energies'].items() \
                                            if index in indiceslist_at_site and index not in exclude and energy is not None] 

            if n_configs is not None: #sort indices from lower to higher energy, and choose the first n_configs
                sorted_items = sorted(indices_and_energies_at_site,key=itemgetter(1))               
                calculations_indices += [x[0] for x in sorted_items[:min(n_configs, len(sorted_items))] ]
            elif threshold is not None:
                e_min = min([x[1] for x in indices_and_energies_at_site])
                calculations_indices += [x[0] for x in indices_and_energies_at_site if x[1] - e_min <= threshold]
    else:
        indices_and_energies = [(index, energy) for index, energy in screening_results['energies'].items() \
                                    if index not in exclude and energy is not None]
        if n_configs is not None: 
            sorted_items = sorted(indices_and_energies,key=itemgetter(1))               
            calculations_indices = [x[0] for x in sorted_items[:min(n_configs, len(sorted_items))] ]
        elif threshold is not None:
            e_min = min([x[1] for x in indices_and_energies])
            calculations_indices = [x[0] for x in indices_and_energies_at_site if x[1] - e_min <= threshold]

    return calculations_indices


def obtain_fullrelax_structures(settings : Settings, calculations_indices : list, REGENERATE : bool = False):
    '''
    Returns a list with the adsorption configurations for the full relax, either reading the final coordinates
    of the screening, or re-generating them from scratch according to settings.in

    Args:
    - settings: Settings object, containing the slab and molecule parameters
    - calculations_indices: indices of the configurations
    - REGENERATE: re-generate the configurations from scratch, starting from isolated slab and molecule. Can be used to change the initial distance
    and repeat the relax from the beginning
    ''' 
    if(REGENERATE):
        all_mol_on_slab_configs_ase, _ = adsorption_configurations(settings=settings)
    else:
        #Read slab and molecule for constraints
        #Slab import from file
        slab = Slab(slab_filename=settings.slab_filename, 
                    layers_threshold=settings.layers_height, 
                    surface_sites_height=settings.surface_height, 
                    fixed_layers_slab=settings.fixed_layers_slab, 
                    fixed_indices_slab=settings.fixed_indices_slab, 
                    fix_slab_xyz=settings.fix_slab_xyz,
                    sort_atoms_by_z=settings.sort_atoms_by_z,
                    translate_slab_from_below_cell_bottom=settings.translate_slab)
        #Molecule import from file
        mol = Molecule(molecule_filename=settings.molecule_filename,
                    atom_index=settings.selected_atom_index,
                    molecule_axis_atoms=settings.molecule_axis_atoms, 
                    axis_vector=settings.axis_vector, 
                    atoms_subset=settings.mol_subset_atoms, 
                    fixed_indices_mol=settings.fixed_indices_mol, 
                    fix_mol_xyz=settings.fix_mol_xyz)        
        
        all_mol_on_slab_configs_ase = []
        for index in calculations_indices:
            atoms = read(OUT_FILE_PATHS['SCREENING'][settings.program].format(index))
            #set constraints from mol and slab if when we are reading the whole pre-relaxed structure
            
            if settings.mol_before_slab:
                c = [constraint for constraint in mol.mol_ase.constraints] + \
                    [FixCartesian(constraint.a + mol.natoms, ~constraint.mask) for constraint in slab.slab_ase.constraints]
                #for some obscure reason FixCartesian wants the mask in the form 1=Fix, 0=Free, but then it negates in internally, so we first need to negate the mask to construct a new object
            else:
                c = [constraint for constraint in slab.slab_ase.constraints] + \
                    [FixCartesian(constraint.a + slab.natoms, ~constraint.mask) for constraint in mol.mol_ase.constraints]              
            
            atoms.set_constraint(c)
            all_mol_on_slab_configs_ase.append(atoms)

    return all_mol_on_slab_configs_ase

 
def generate(SAVEFIG=False):
    '''
    Generates adsorption configurations and writes the inputs, using the settings for the final relaxations.
    Useful to use the program just as a generator of configurations

    Args:
    - SAVEFIG: save an image of the adsorption sites and of the molecular rotations when generating the configurations
    '''
    
    settings=Settings(read_energies=False)

    all_mol_on_slab_configs_ase, full_labels = adsorption_configurations(settings, SAVEFIG)

    write_labels_csvfile(full_labels, labels_filename=labels_filename)

    indices_list = np.arange(len(all_mol_on_slab_configs_ase))

    write_inputs(settings, 
                 all_mol_on_slab_configs_ase,
                 indices_list,
                 calc_type='SCREENING',
                 OVERRIDE_SETTINGS=False,
                 INTERACTIVE=True)

  
def launch_screening(SAVEFIG : bool = False):
    '''
    Generates adsorption configurations, writes inputs and launches calculations for the preliminary screening.

    Args:
    - SAVEFIG: save an image of the adsorption sites and of the molecular rotations when generating the configurations
    '''
    
    settings=Settings()

    all_mol_on_slab_configs_ase, full_labels = adsorption_configurations(settings, SAVEFIG, VERBOSE=True)

    write_labels_csvfile(full_labels, labels_filename=labels_filename)

    written_indices = write_inputs(settings, 
                                    all_mol_on_slab_configs_ase,
                                    calc_indices=np.arange(len(all_mol_on_slab_configs_ase)),
                                    calc_type='SCREENING', 
                                    INTERACTIVE=True)

    launch_jobs(program=settings.program,
                calc_type='SCREENING',
                jobscript=settings.jobscript,
                sbatch_command=settings.sbatch_command,
                indices_list=written_indices)    

 
def final_relax(n_configs: int = None, threshold : float = None, exclude : list= None, required_indices : list = None, REGENERATE=False, BY_SITE = False):
    '''
    Reads/generates adsorption configurations, writes inputs and launches calculations for the final relax.

    Args:
    - n_configs: nubmer of configurations to be relaxed, starting from the one with lowest energy
    - threshold: energy threshold (in eV) from the NOT EXCLUDED lowest energy configuration. The configuration with E - Emin < threshold will be selected
    - exclude: indices of the configurations to be excluded
    - required_indices: user-specified indices, instead of identifying them according to energys
    - REGENERATE: re-generate the configurations from scratch, starting from isolated slab and molecule. Can be used to change the initial distance
    and repeat the relax from the beginning
    - BY_SITE: do the configuration identification separately for each site. One or more configuration for each site will be produced
    '''

    #Check input parameters and set n_configs accordingly
    if n_configs is None and threshold is None: n_configs = N_relax_default if not BY_SITE else 1

    settings=Settings()

    if required_indices: calculations_indices = required_indices
    else: 
        calculations_indices = obtain_fullrelax_indices(settings=settings, 
                                                        n_configs=n_configs, 
                                                        threshold=threshold, 
                                                        exclude=exclude, 
                                                        BY_SITE=BY_SITE)

    all_mol_on_slab_configs_ase = obtain_fullrelax_structures(settings, calculations_indices, REGENERATE)


    write_inputs(settings, 
                all_mol_on_slab_configs_ase,
                calculations_indices,
                calc_type='RELAX', 
                INTERACTIVE=True)

    launch_jobs(program=settings.program,
                calc_type='RELAX',
                jobscript=settings.jobscript,
                sbatch_command=settings.sbatch_command,
                indices_list=calculations_indices)   
