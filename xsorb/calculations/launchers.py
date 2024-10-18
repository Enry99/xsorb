"""
Created on Tue 28 Feb 2023
@author: Enrico Pedretti

Main functions to launch the various types of calculations.

"""

import sys

import numpy as np

from xsorb.structures import AdsorptionStructuresGenerator
from xsorb.io.settings import Settings
from xsorb.io.utils import ase_custom_read as read
from xsorb.io.inputs import write_inputs, write_slab_mol_inputs
from xsorb.io.jobs import launch_jobs
from xsorb.io.database import Database
from xsorb.io.utils import continue_even_if_not_all_completed_question
from xsorb.calculations.selection import obtain_calc_indices, get_adsorption_structures

N_RELAX_DEFAULT = 5


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
        calc_ids = obtain_calc_indices(calc_type='ml_opt')
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
                jobscript=settings.input.jobscript_ml_path,
                sbatch_command=settings.input.submit_command_ml,
                systems=written_systems,
                jobname_prefix=settings.input.jobname_prefix)



def launch_final_relax(*,
                       n_configs: int | None = None,
                       threshold : float | None = None,
                       calc_ids : list[int] | None = None,
                       excluded_calc_ids : list[int] | None = None,
                       take_from : str = 'screening',
                       relax_from_initial : bool = False,
                       by_site : bool = False,
                       by_mol_idx : bool = False,
                       separate_chem_phys : bool = False):
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
    - by_mol_atom: do the configuration identification separately for each molecule ref. atom.
    - separate_chem_phys: do the configuration identification separately
        for physisorption and chemisorption
    '''

    #Initial setup of parameters

    if take_from not in ('screening', 'ml_opt'):
        raise ValueError('Invalid value for take_from. Must be "screening" or "ml_opt".')

    #check that only one between n_configs, threshold, required_calc_ids is specified,
    if np.sum([n_configs is not None, threshold is not None, calc_ids is not None]) > 1:
        raise RuntimeError('Only one between n_configs, threshold, '\
                           'required_calc_ids can be specified.')
    if n_configs is None and threshold is None and calc_ids is None:
        #none specified, use n_configs method, with default values
        if by_site:
            n_configs = 1
        else:
            n_configs = N_RELAX_DEFAULT


    #Retrieve the structures
    settings=Settings()

    #this check also updates the db
    if not Database.all_completed(calc_type=take_from) and \
        not continue_even_if_not_all_completed_question():
        print('Quitting.')
        sys.exit(0)

    if not calc_ids:
        #retrieve the indices of the configurations to be relaxed
        calc_ids = obtain_calc_indices(calc_type=take_from,
                                       n_configs=n_configs,
                                       threshold=threshold,
                                       excluded_calc_ids=excluded_calc_ids,
                                       by_site=by_site,
                                       by_mol_atom=by_mol_idx,
                                       separate_chem_phys=separate_chem_phys)
    else:
        #use the user-specified indices, exclude unwanted calculations
        if excluded_calc_ids:
            calc_ids = [calc_id for calc_id in calc_ids if calc_id not in excluded_calc_ids]

    get_structures_from = 'structures' if relax_from_initial else take_from
    adsorption_structures = get_adsorption_structures(get_structures_from, calc_ids)

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


def launch_isolated_slab_and_molecule(ml : bool, mol_only : bool = False):
    '''
    Launch the calculations for the isolated slab and molecule.

    Args:
    - ml: for machine learning calculations
    '''

    settings = Settings()

    slab = read(settings.input.slab_filename) if not mol_only else None
    mol = read(settings.input.molecule_filename)

    written_systems = write_slab_mol_inputs(slab=slab, molecule=mol, settings=settings, ml=ml)

    launch_jobs(program=settings.program,
                calc_type='isolated',
                jobscript=settings.input.jobscript_path if not ml \
                    else settings.input.jobscript_ml_path,
                sbatch_command=settings.input.submit_command if not ml \
                    else settings.input.submit_command_ml,
                systems=written_systems,
                jobname_prefix=settings.input.jobname_prefix)
