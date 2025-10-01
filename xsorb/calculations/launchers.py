#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Enrico Pedretti

"""
Created on Tue 28 Feb 2023
@author: Enrico Pedretti

Main functions to launch the various types of calculations.

"""

from __future__ import annotations
import sys
import logging

import numpy as np

from xsorb.structures.generation import AdsorptionStructuresGenerator
from xsorb.settings import Settings
from xsorb.ase_custom.io import ase_custom_read as read
from xsorb.io.inputs import write_inputs, write_slab_mol_inputs
from xsorb.io.jobs import launch_jobs
from xsorb.io.cleanup import fresh_start
from xsorb.io.database import Database
from xsorb.io.utils import yes_no_question
from xsorb.calculations.selection import obtain_calc_indices, get_adsorption_structures


def generate(save_image : bool = False):
    '''
    Generates adsorption configurations and writes the inputs,
    using the settings for the final relaxations.
    Useful to inspect the generated configurations before launching the calculations.

    Args:
    - save_image: save an image of the adsorption sites
    and of the molecular rotations after writing the files.
    '''

    fresh_start()  # ask the user if they want to start a fresh run

    settings=Settings()

    slab = read(settings.input.slab_filename)
    mol = read(settings.input.molecule_filename)

    gen = AdsorptionStructuresGenerator(slab, mol, settings, verbose=True)
    adsorption_structures = gen.generate_adsorption_structures(save_image=save_image)

    write_inputs(adsorption_structures=adsorption_structures, settings=settings)



def launch_screening(from_ml_opt : bool = False, save_image : bool = False,):
    '''
    Generates adsorption configurations, writes inputs
    and launches calculations for the preliminary screening.

    Args:
    - from_ml_opt: use the final configuration from machine learning pre-optimization
        as starting point instead of generating them from scratch
    - save_image: save an image of the adsorption sites
    and of the molecular rotations when generating the configurations
    '''

    if not from_ml_opt:
        fresh_start()

    settings=Settings(read_energies_dft=True) #need the energies to store them into the db metadata

    if from_ml_opt:
    #this check also updates the db
        if not Database.all_completed(calc_type='mlopt') and \
            not yes_no_question('Not all calculations are completed. '\
                                'Continue anyway with those that are present?'):
            logging.info('Quitting.')
            sys.exit(0)
        calc_ids = obtain_calc_indices(calc_type='mlopt')
        adsorption_structures = get_adsorption_structures(calc_ids=calc_ids,
                                                          get_structures_from='mlopt')
    else:
        calc_ids = None # generate the ids when inserting the structures into the database
        slab = read(settings.input.slab_filename)
        mol = read(settings.input.molecule_filename)
        gen = AdsorptionStructuresGenerator(slab, mol, settings, verbose=True)
        adsorption_structures = gen.generate_adsorption_structures(save_image=save_image)


    written_systems = write_inputs(adsorption_structures=adsorption_structures,
                                   settings=settings,
                                   calc_type='screening',
                                   calc_ids=calc_ids)

    launch_jobs(program=settings.dft.program,
                calc_type='screening',
                jobscript=settings.input.jobscript_path,
                scheduler_name=settings.input.scheduler,
                systems_calcinfos=written_systems,
                jobname_prefix=settings.input.jobname_prefix)



def launch_ml_opt(save_image : bool = False,):
    '''
    Generates adsorption configurations, writes inputs
    and launches calculations for the machine learning optimization.
    Args:
    - save_image: save an image of the adsorption sites
    and of the molecular rotations when generating the configurations
    '''

    fresh_start()  # ask the user if they want to start a fresh run

    settings=Settings(read_energies_ml=True) #need the energies to store them into the db metadata

    slab = read(settings.input.slab_filename)
    mol = read(settings.input.molecule_filename)

    gen = AdsorptionStructuresGenerator(slab, mol, settings, verbose=True)
    adsorption_structures = gen.generate_adsorption_structures(save_image=save_image)

    written_systems = write_inputs(adsorption_structures=adsorption_structures,
                                   settings=settings,
                                   calc_type='mlopt')

    if settings.input.jobscript_ml_path is None:
        raise ValueError('jobscript_ml_path is not defined in the settings file. '\
                         'Please define it to launch the machine learning optimization.')

    launch_jobs(program='ml',
                calc_type='mlopt',
                jobscript=settings.input.jobscript_ml_path,
                scheduler_name=settings.input.scheduler,
                systems_calcinfos=written_systems,
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
    - calc_ids: user-specified indices, instead of identifying them according to energy
    - excluded_calc_ids: indices of the configurations to be excluded
    - take_from: type of calculation for the selection. Can be 'screening', 'mlopt'
    - relax_from_initial: use the initial configuration as starting point for the relaxation
    - by_site: do the configuration identification separately for each site.
        One or more configuration for each site will be produced
    - by_mol_atom: do the configuration identification separately for each molecule ref. atom.
    - separate_chem_phys: do the configuration identification separately
        for physisorption and chemisorption
    '''

    #Initial setup of parameters

    if take_from not in ('screening', 'mlopt'):
        raise ValueError('Invalid value for take_from. Must be "screening" or "mlopt".')

    #check that only one between n_configs, threshold, required_calc_ids is specified,
    if np.sum([n_configs is not None, threshold is not None, calc_ids is not None]) > 1:
        raise RuntimeError('Only one between n_configs, threshold, '\
                           'required_calc_ids can be specified.')
    if n_configs is None and threshold is None and calc_ids is None:
        raise RuntimeError(
            'At least one of n_configs, threshold, required_calc_ids must be specified.')


    #Retrieve the structures
    settings=Settings(read_energies_dft=True) #need the energies to store them into the db metadata

    #this check also updates the db
    if not Database.all_completed(calc_type=take_from) and \
        not yes_no_question('Not all calculations are completed. '\
                            'Continue anyway with those that are present?'):
        logging.info('Quitting.')
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

    if len(calc_ids) == 0:
        logging.info('No configurations to be relaxed. Quitting.')
        sys.exit(0)

    get_structures_from = 'structures' if relax_from_initial else take_from
    adsorption_structures = get_adsorption_structures(get_structures_from, calc_ids)

    written_systems = write_inputs(adsorption_structures=adsorption_structures,
                                   settings=settings,
                                   calc_type='relax',
                                   calc_ids=calc_ids)

    launch_jobs(program=settings.dft.program,
                calc_type='relax',
                jobscript=settings.input.jobscript_path,
                scheduler_name=settings.input.scheduler,
                systems_calcinfos=written_systems,
                jobname_prefix=settings.input.jobname_prefix)



def launch_isolated_slab_and_molecule(*,
                                      ml : bool,
                                      launch_slab : bool = True,
                                      launch_mol: bool = True,
                                      samecell : bool = False,
                                      use_constraints : bool = False):
    '''
    Launch the calculations for the isolated slab and molecule.

    Args:
    - ml: for machine learning calculations
    - launch_slab: launch the calculations for the slab
    - launch_mol: launch the calculations for the molecule
    - samecell: use the same slab cell also for the molecule (to remove coverage effects)
    - use_constraints: include the constraints defined in settings file
    '''

    settings = Settings()

    slab = read(settings.input.slab_filename)
    mol = read(settings.input.molecule_filename)

    if slab.cell is None:
        raise RuntimeError('The slab cell is not defined.')

    if samecell:
        mol.cell = slab.cell
    elif not mol.cell:
        mol.center(vacuum=5.0)

    slab.pbc = True
    mol.pbc = True

    if use_constraints:
        from xsorb.structures.slab import Slab
        from xsorb.structures.molecule import Molecule

        slab = Slab(slab=slab,
                    layers_threshold=settings.structure.constraints.layers_height,
                    fixed_layers_slab=settings.structure.constraints.fixed_layers_slab,
                    fixed_indices_slab=settings.structure.constraints.fixed_indices_slab,
                    fixed_thickness_slab=settings.structure.constraints.fixed_thickness_slab,
                    fix_slab_xyz=settings.structure.constraints.fix_slab_xyz
                    ).slab_ase

        mol = Molecule(mol=mol,
                    atom_indexes=settings.structure.molecule.selected_atom_indexes,
                    molecule_axis_mode=settings.structure.molecule.molecule_axis.mode,
                    molecule_axis_values=settings.structure.molecule.molecule_axis.values,
                    fixed_indices_mol=settings.structure.constraints.fixed_indices_mol,
                    fix_mol_xyz=settings.structure.constraints.fix_mol_xyz
                    ).mol_ase


    written_systems = write_slab_mol_inputs(slab=slab if launch_slab else None,
                                            molecule=mol if launch_mol else None,
                                            settings=settings,
                                            ml=ml,
                                            force_gamma=not samecell)

    if ml:
        if settings.input.jobscript_ml_path is None:
            raise ValueError('jobscript_ml_path is not defined in the settings file. '\
                'Please define it to launch the machine learning optimization.')
        program = 'ml'
        jobscript = settings.input.jobscript_ml_path
    else:
        program = settings.dft.program
        jobscript = settings.input.jobscript_path

    launch_jobs(program=program,
                calc_type='isolated',
                jobscript=jobscript,
                scheduler_name=settings.input.scheduler,
                systems_calcinfos=written_systems,
                jobname_prefix=settings.input.jobname_prefix)
