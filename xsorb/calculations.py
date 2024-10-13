"""
Created on Tue 28 Feb 2023
@author: Enrico Pedretti

Main functions to generate the adsorption configurations and set up the calculations

"""

import os
import sys
from operator import itemgetter

import numpy as np
import pandas as pd
from ase.io import read, write
from ase.constraints import FixCartesian


from xsorb.settings import Settings
from xsorb.dft_codes.definitions import OUT_FILE_PATHS, IN_FILE_PATHS
from xsorb.common_definitions import *
from xsorb.structures import AdsorptionStructuresGenerator, AdsorptionStructure
from xsorb.io.inputs import write_inputs
from xsorb.io.launch import launch_jobs


def regenerate_missing_sitelabels():
    '''
    Re-generate the site_labels.csv file if it was accidentally deleted, using the info from settings.in
    '''

    settings=Settings()

    _, full_labels, _, _ = adsorption_configurations(settings, save_image=False, VERBOSE=True)

    write_labels_csvfile(full_labels, labels_filename=labels_filename)




def write_slab_mol_ml(slab, mol):
    os.makedirs(preopt_outdir, exist_ok=True)
    os.makedirs(preopt_outdir+'/slab', exist_ok=True)
    os.makedirs(preopt_outdir+'/mol', exist_ok=True)

    write(preopt_outdir+'/slab/slab.xyz', slab)
    write(preopt_outdir+'/mol/mol.xyz', mol)


def obtain_fullrelax_indices(settings : Settings,
                             n_configs: int = None,
                             threshold : float = None,
                             exclude : list= None,
                             BY_SITE = False,
                             from_preopt : bool = False):
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
    if from_preopt:
        screening_results = get_calculations_results_ml()
    else:
        screening_results = get_calculations_results(program=settings.program, calc_type='SCREENING')
    print('Screening energies collected.')


    #Ask to quit if some screening calculations are missing
    if False in screening_results['relax_completed'].values() \
        or len(screening_results['relax_completed'].values()) < len(_get_configurations_numbers()):

        print('Not all screening calculations have provided a final energy." \
                "Those for which no output exists will be excluded, while for those not completed the last coordinates will be used.')
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
            calculations_indices = [x[0] for x in indices_and_energies if x[1] - e_min <= threshold]

    return calculations_indices


def obtain_fullrelax_structures(settings : Settings, calculations_indices : list, REGENERATE : bool = False, from_preopt : bool = False):
    '''
    Returns a list with the adsorption configurations for the full relax, either reading the final coordinates
    of the screening, or re-generating them from scratch according to settings.in

    Args:
    - settings: Settings object, containing the slab and molecule parameters
    - calculations_indices: indices of the configurations
    - REGENERATE: re-generate the configurations from scratch, starting from isolated slab and molecule. Can be used to change the initial distance
    and repeat the relax from the beginning
    '''

    all_mol_on_slab_configs_ase, _, _, _ = adsorption_configurations(settings=settings)

    if REGENERATE:
        return all_mol_on_slab_configs_ase
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

        new_all_mol_on_slab_configs_ase = []
        for index in calculations_indices:
            atoms = read(OUT_FILE_PATHS['SCREENING' if not from_preopt else 'PREOPT'][settings.program if not from_preopt else 'ML'].format(index))
            #set constraints from mol and slab if when we are reading the whole pre-relaxed structure

            if settings.mol_before_slab:
                c = [constraint for constraint in mol.mol_ase.constraints] + \
                    [FixCartesian(constraint.a + mol.natoms, ~constraint.mask) for constraint in slab.slab_ase.constraints]
                #for some obscure reason FixCartesian wants the mask in the form 1=Fix, 0=Free, but then it negates in internally, so we first need to negate the mask to construct a new object
            else:
                c = [constraint for constraint in slab.slab_ase.constraints] + \
                    [FixCartesian(constraint.a + slab.natoms, ~constraint.mask) for constraint in mol.mol_ase.constraints]

            origatoms = all_mol_on_slab_configs_ase[index]
            origatoms.set_constraint(c)
            origatoms.positions = atoms.positions
            new_all_mol_on_slab_configs_ase.append(origatoms)

        return new_all_mol_on_slab_configs_ase


def generate(save_image : bool = False):
    '''
    Generates adsorption configurations and writes the inputs,
    using the settings for the final relaxations.
    Useful to inspect the generated configurations before launching the calculations.

    Args:
    - save_image: save an image of the adsorption sites
    and of the molecular rotations after writing the files.
    '''

    settings=Settings()

    gen = AdsorptionStructuresGenerator(settings, verbose=True)
    adsorption_structures = gen.generate_adsorption_structures(write_sites=False,
                                                               save_image=save_image)

    write_inputs(adsorption_structures=adsorption_structures, settings=settings)


def preopt_ml(save_image : bool = False):

    settings=Settings()


    all_mol_on_slab_configs_ase, full_labels, slab, mol = adsorption_configurations(settings, save_image, VERBOSE=True)

    write_labels_csvfile(full_labels, labels_filename=labels_filename)

    write_slab_mol_ml(slab=slab, mol=mol)

    written_indices = write_inputs_ml(settings,
                                    all_mol_on_slab_configs_ase,
                                    calc_indices=np.arange(len(all_mol_on_slab_configs_ase))
                                    )

    #slab and molecule
    launch_jobs_ml(jobscript_ml=settings.jobscript_ml,
                   sbatch_command=settings.sbatch_command_ml,
                    explicit_labels=['slab', 'mol'],
                    fix_bondlengths=settings.fix_bondlengths_preopt,
                    fix_slab=settings.fix_slab_preopt,
                    slab_indices=[0 + settings.mol_before_slab * len(mol), len(slab) + settings.mol_before_slab * len(mol)],
                    jobname_prefix=settings.jobname_prefix)

    #adsorption structures
    launch_jobs_ml(jobscript_ml=settings.jobscript_ml,
                sbatch_command=settings.sbatch_command_ml,
                indices_list=written_indices,
                fix_bondlengths=settings.fix_bondlengths_preopt,
                fix_slab=settings.fix_slab_preopt,
                slab_indices=[0 + settings.mol_before_slab * len(mol), len(slab) + settings.mol_before_slab * len(mol)],
                jobname_prefix=settings.jobname_prefix)


def get_preopt_structures():
    '''
    Reads the preopt structures from the output files of the preopt calculations
    '''
    #settings=Settings()

    config_indices = _get_configurations_numbers()

    all_mol_on_slab_configs_ase = []
    for index in config_indices:
        atoms = read(OUT_FILE_PATHS['PREOPT']['ML'].format(index))
        all_mol_on_slab_configs_ase.append(atoms)

    return all_mol_on_slab_configs_ase


def launch_screening(from_preopt : bool = False, save_image : bool = False,):
    '''
    Generates adsorption configurations, writes inputs
    and launches calculations for the preliminary screening.

    Args:
    - from_preopt: use the final configuration from machine learning pre-optimization
        as starting point instead of generating them from scratch
    - save_image: save an image of the adsorption sites
    and of the molecular rotations when generating the configurations
    '''

    settings=Settings()

    if from_preopt:
        adsorption_structures = get_preopt_structures()
    else:
        gen = AdsorptionStructuresGenerator(settings, verbose=True)
        adsorption_structures = gen.generate_adsorption_structures(write_sites=True,
                                                                    save_image=save_image)

    written_systems = write_inputs(adsorption_structures=adsorption_structures, settings=settings)


    launch_jobs(program=settings.program,
                calc_type='screening',
                jobscript=settings.input.jobscript_path,
                sbatch_command=settings.input.submit_command,
                systems=written_systems,
                jobname_prefix=settings.input.jobname_prefix)


def final_relax(n_configs: int = None, threshold : float = None, exclude : list= None, required_indices : list = None, from_preopt : bool = False, REGENERATE=False, BY_SITE = False):
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
                                                        BY_SITE=BY_SITE,
                                                        from_preopt=from_preopt)

    all_mol_on_slab_configs_ase = obtain_fullrelax_structures(settings, calculations_indices, REGENERATE, from_preopt=from_preopt)


    write_inputs(settings,
                all_mol_on_slab_configs_ase,
                calculations_indices,
                calc_type='RELAX',
                INTERACTIVE=True)

    launch_jobs(program=settings.program,
                calc_type='RELAX',
                jobscript=settings.jobscript,
                sbatch_command=settings.sbatch_command,
                indices_list=calculations_indices,
                jobname_prefix=settings.jobname_prefix)
