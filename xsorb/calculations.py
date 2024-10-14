"""
Created on Tue 28 Feb 2023
@author: Enrico Pedretti

Main functions to generate the adsorption configurations and set up the calculations

"""

import sys
from operator import itemgetter

import numpy as np
from ase import Atoms


from xsorb.settings import Settings
from xsorb.common_definitions import *
from xsorb.structures import AdsorptionStructuresGenerator, AdsorptionStructure
from xsorb.io.inputs import write_inputs
from xsorb.io.launch import launch_jobs
from xsorb.io.database import Database
from xsorb.io.utils import continue_even_if_not_all_completed_question





def obtain_calc_indices(calc_type : str,
                        n_configs: int = None,
                        threshold : float = None,
                        excluded_calc_ids : list= None,
                        by_site : bool = False) -> list[int]:
    '''
    Returns a list with the indices of the configurations to be relaxed, according to the specified criteria.

    Args:
    - calc_type: type of calculation to get the indices from. Can be 'screening', 'ml_opt'
    - n_configs: number of configurations to be relaxed, starting from the one with lowest energy
    - threshold: energy threshold (in eV) from the NOT EXCLUDED lowest energy configuration.
        The configuration with E - Emin < threshold will be selected
    - excluded_calc_ids: indices of the configurations to be excluded
    - by_site: do the configuration identification separately for each site.

    Returns:
    - calculations_indices: list of indices of the configurations to be relaxed
    '''


    print(f'Collecting results from {calc_type}...')
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


def get_adsorption_structures(calc_ids : list[int] | None = None,
                              get_structures_from : bool = False) -> list[AdsorptionStructure]:
    '''
    Returns the a list of AdsorptionStructure objects, but with the atoms
    substituted with the ones from the previous calculation.
    If from ml_opt, the the original constraints are set, important when
    the ml_opt was performed by fixing the slab.

    Args:
    - calc_ids: list of indices of the configurations to be retrieved. If None, all are retrieved
    - get_structures_from: type of calculation to get the structures from.
        Can be 'screening', 'ml_opt', 'structures'

    Returns:
    - adsorption_structures: list of AdsorptionStructure objects
    '''

    #get structures from database
    rows = Database.get_calculations(calc_type=get_structures_from, calc_ids=calc_ids)

    if get_structures_from == 'ml_opt':
        rows_original = Database.get_calculations(calc_type='structures', calc_ids=calc_ids)
        constraints = [row.constraints for row in rows_original]

    # prepare the structures by substituting the atoms with the one from the
    # previous calculation
    adsorption_structures : list [AdsorptionStructure] = []
    for row in rows:
        ads_struct : AdsorptionStructure = row.data.adsorption_structure
        atoms : Atoms = row.toatoms()

        if get_structures_from == 'ml_opt':
            atoms.set_constraint(constraints.pop(0))

        ads_struct.atoms = atoms
        adsorption_structures.append(ads_struct)

    return adsorption_structures


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


def launch_isolated_slab_and_molecule(program_type : str):
    '''
    Launch the calculations for the isolated slab and molecule.

    Args:
    - program_type: 'dft' or 'ml'
    '''
    pass


def launch_ml_opt(save_image : bool = False,):
    '''
    Generates adsorption configurations, writes inputs
    and launches calculations for the machine learning optimization.
    Args:
    - save_image: save an image of the adsorption sites
    and of the molecular rotations when generating the configurations
    '''

    settings=Settings()

    gen = AdsorptionStructuresGenerator(settings, verbose=True)
    adsorption_structures = gen.generate_adsorption_structures(write_sites=True,
                                                                save_image=save_image)

    written_systems = write_inputs(adsorption_structures=adsorption_structures,
                                   settings=settings,
                                   calc_type='ml_opt')


    launch_jobs(program=settings.program,
                calc_type='ml_opt',
                jobscript=settings.input.jobscript_path,
                sbatch_command=settings.input.submit_command,
                systems=written_systems,
                jobname_prefix=settings.input.jobname_prefix)


def launch_screening(from_ml_opt : bool = False, save_image : bool = False,):
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

    if from_ml_opt:
        calc_ids = #TODO:WRITE
        adsorption_structures = get_adsorption_structures(calc_ids=calc_ids,
                                                          get_structures_from='ml_opt')
    else:
        gen = AdsorptionStructuresGenerator(settings, verbose=True)
        adsorption_structures = gen.generate_adsorption_structures(write_sites=True,
                                                                    save_image=save_image)
        calc_ids = None


    written_systems = write_inputs(adsorption_structures=adsorption_structures,
                                   settings=settings,
                                   calc_type='screening',
                                   calc_ids=calc_ids)


    launch_jobs(program=settings.program,
                calc_type='screening',
                jobscript=settings.input.jobscript_path,
                sbatch_command=settings.input.submit_command,
                systems=written_systems,
                jobname_prefix=settings.input.jobname_prefix)


def launch_final_relax(*,
                       n_configs: int = None,
                       threshold : float = None,
                       calc_ids : list[int] = None,
                       excluded_calc_ids : list[int] = None,
                       take_from : str = 'screening',
                       relax_from_initial : bool = False,
                       by_site : bool = False):
    '''
    Reads/generates adsorption configurations, writes inputs and launches the
    calculations for the final relax. If neither n_configs, threshold, required_calc_ids
    is specified, the n_configs mode is used with default values.

    Args:
    - n_configs: nubmer of configurations to be relaxed, starting from the one with lowest energy
    - threshold: energy threshold (in eV) from the NOT EXCLUDED lowest energy configuration.
        The configuration with E - Emin < threshold will be selected
    - required_calc_ids: user-specified indices, instead of identifying them according to energy
    - excluded_calc_ids: indices of the configurations to be excluded
    - take_from: type of calculation for the selection. Can be 'screening', 'ml_opt'
    - relax_from_initial: use the initial configuration as starting point for the relaxation
    - by_site: do the configuration identification separately for each site.
        One or more configuration for each site will be produced
    '''

    #Initial setup of parameters

    if take_from not in ('screening', 'ml_opt'):
        raise ValueError('Invalid value for take_from. Must be "screening" or "ml_opt".')

    #check that only one between n_configs, threshold, required_calc_ids is specified,
    if np.sum([n_configs is not None, threshold is not None, calc_ids is not None]) > 1:
        raise RuntimeError('Only one between n_configs, threshold, '\
                           'required_calc_ids can be specified.')
    elif n_configs is None and threshold is None and calc_ids is None:
        #none specified, use n_configs method, with default values
        if by_site:
            n_configs = 1
        else:
            n_configs = N_relax_default


    #Retrieve the structures
    settings=Settings()

    if not Database.all_completed(calc_type=take_from) and \
        not continue_even_if_not_all_completed_question():
        print('Quitting.')
        sys.exit(0)

    if not calc_ids:
        #retrieve the indices of the configurations to be relaxed

        calc_ids = [] #TODO:WRITE
    else:
        #use the user-specified indices, exclude unwanted calculations
        if excluded_calc_ids:
            calc_ids = [calc_id for calc_id in calc_ids if calc_id not in excluded_calc_ids]

    get_structures_from = 'structures' if relax_from_initial else take_from
    adsorption_structures = get_adsorption_structures(calc_ids, get_structures_from)


    written_systems = write_inputs(adsorption_structures=adsorption_structures,
                                   settings=settings,
                                   calc_type='relax',
                                   calc_ids=calc_ids)


    launch_jobs(program=settings.program,
                calc_type='relax',
                jobscript=settings.input.jobscript_path,
                sbatch_command=settings.input.submit_command,
                systems=written_systems,
                jobname_prefix=settings.input.jobname_prefix)
