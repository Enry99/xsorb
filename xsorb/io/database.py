#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Enrico Pedretti

'''
Module to handle the ase database

Four databases are used:
- structures.json: Database with the initial structures
- screening.json: Contains the results of the screening calculations
- relaxations.json: Contains the results of the relaxations
- mlopt.json: Contains the results of the machine learning optimization


The calculation databases have the following columns:

--- Default columns from ase.db ---
- id (int): unique id of the row
- atoms (Atoms): Atoms object with the structure
- energy (float): total energy of the structure (optional)

--- Extra columns, defined in xsorb.adsorptiondata.AdsorptionCalculation.db_keys() ---
AdsorptionStructure data:
- "site", "site_info", "xrot", "yrot", "zrot", "mol_atom", "conform_id", "coords", "initial_dz"

CalculationInfo data:
- calc_id (str): unique id of the calculation (consistent with numeration in structures.json)
- in_file_path (str): path to the input file
- out_file_path (str): path to the output file
- log_file_path (str): path to the log file

CalculationResults data:
- adsorption_energy (float): adsorption energy of the structure
- status (str): status of the calculation (incomplete, completed, scf_nonconverged)
- bonds (str): string with the bonds between the molecule and the slab
- final_dz (float): final vertical distance between the reference atom of the molecule
    and the adsorption site

--- Extra columns to handle job submission ---
- job_id (str): id of the job in the scheduler


--- data (dict) from ase, used to store more complex data structures ---
- data contains a dictionary with xsorb.adsorptiondata.AdsorptionCalculation,
    which can be converted back to an object


Each calculation database also has the following metadata:
- program (str): name of the program used for the calculations
- mult (float): multiplicative factor for the covalent radii to determine bonding
    The mult factor is updated from settings only when refreshing the database
- e_slab (float): total energy of the isolated slab
- e_mol (float or list of float): total energy of the isolated molecule
    (if multiple conformers are used, this is a list)
- adsorption_sites (list): list of AdsorptionSiteCrystal or AdsorptionSiteAmorphous objects

'''

from __future__ import annotations
from pathlib import Path
import logging

import pandas as pd
import ase.db

import xsorb.calculations.results
from xsorb.io.filenames import STRUCTURES_DB_NAME, CALC_DB_NAMES
from xsorb.ase_custom.atoms import AtomsCustom
from xsorb.adsorptiondata.adsorptionstructure import (
    AdsorptionSiteCrystal, AdsorptionSiteAmorphous, SurroundingSite)
from xsorb.adsorptiondata import AdsorptionCalculation
from xsorb.adsorptiondata import AdsorptionStructure
from xsorb.settings.settings import Energies


class Database:
    '''
    Collection of static functions to read/write generated structures and calculations
    results to ase db
    '''


    @staticmethod
    def add_structures(adsorption_structures : list [AdsorptionStructure],
                       write_csv : bool = True) -> list[int]:
        '''
        Adds the generated adsorption structures to the structures database,
        returning the corresponding calc_ids

        Args:
        - adsorption_structures: list of AdsorptionStructure objects
        - write_csv: write the info to a csv file

        Returns:
        - list[int] with the row ids of the structures,
            which will become the calc_ids for the calculations
        '''

        # Write the adsorption structures to the database
        # excluding those that are already present
        calc_ids : list[int] = []
        with ase.db.connect(STRUCTURES_DB_NAME) as db:
            try:
                metadata = db.metadata.copy() # workaround for ase.db bug
            except: #pylint: disable=bare-except
                metadata = {}
            adsorption_sites = Database.get_adsorption_sites() # get existing sites
            for ads_struct in adsorption_structures:
                already_present = False
                for row in db.select(include_data=False):
                    if ads_struct.atoms == AtomsCustom(row.toatoms()):
                        already_present = True
                        calc_ids.append(row.id)
                        break
                if not already_present:
                    calc_id = db.write(ads_struct.atoms,
                                       data={'AdsorptionStructure': ads_struct},
                                       **ads_struct.db_keys())
                    calc_ids.append(calc_id)

                # add the adsorption sites if not already present
                # do not write SurroundingSites, as the they are
                # already contained in amorphous main sites
                if ads_struct.adsite not in adsorption_sites and \
                    not isinstance(ads_struct.adsite, SurroundingSite):
                    adsorption_sites.append(ads_struct.adsite)

            # Update the metadata with the adsorption sites
            metadata['adsorption_sites'] = [site for site in adsorption_sites]

            db.metadata = metadata # workaround for ase.db bug

        if write_csv:
            Database.write_csvfile(include_results=False)

        return calc_ids


    @staticmethod
    def add_calculations(systems : list [AdsorptionCalculation],
                         program : str,
                         mult : float,
                         store_full_trajectories : bool,
                         energies : Energies,
                         calc_type : str) -> None:
        '''
        Write new calculations to corresponding database

        Args:
        - systems: list of AdsorptionCalculation objects, each containing
            adsorption_structure, (calc_info, calc_results)
        - program: string with the name of the program used for the calculations:
            'vasp', 'espresso', 'ml'
        - mult: float with the multiplicative factor for the covalent radii to
            determine bonding
        - energies: Energies object with slab and mol energies
        - calc_type: screening, relax or mlopt
        '''

        # Write the adsorption structures to the corresponding database
        with ase.db.connect(CALC_DB_NAMES[calc_type]) as db:

            for system in systems:
                try:
                    #if the calculation is already present in the database, remove it
                    del db[db.get(f'calc_id={int(system.calc_info.calc_id)}').id]
                except KeyError:
                    #if the calculation is not present, do nothing
                    pass

                #write the new calculation in any case
                ads_struct : AdsorptionStructure = system.adsorption_structure
                db.write(ads_struct.atoms,
                        data={'AdsorptionCalculation': system},
                        **system.db_keys())

            db.metadata = {'program': program,
                           'mult': mult,
                           'store_full_trajectories': store_full_trajectories,
                           'e_slab': energies.E_slab_ml if program == 'ml' else energies.E_slab,
                           'e_mol': energies.E_mol_ml if program == 'ml' else energies.E_mol}


    @staticmethod
    def get_updated_energies(metadata : dict,
                             energies : Energies,
                             program : str,
                             calc_type : str) -> Energies:
        '''
        Get the slab and molecule energies, either from the provided Energies object,
        or from the database metadata if not provided.

        Parameters
        ----------
        - metadata: dict with the database metadata
        - energies: Energies object with slab and mol energies, or None
        - program: string with the name of the program used for the calculations:
            'vasp', 'espresso', 'ml'
        - calc_type: string with the type of calculation: 'screening', 'relax', 'mlopt'

        Returns
        ----------
        - Energies object with the slab and mol energies
        '''

        # hack workaround for backward compatibility
        if 'total_e_slab_mol' in metadata:
            if program == 'ml':
                return Energies(E_slab_ml=metadata['total_e_slab_mol'], E_mol_ml=0.0)
            else:
                return Energies(E_slab=metadata['total_e_slab_mol'], E_mol=0.0)

        if energies is not None:
            # store the energies in the metadata and return the original object
            if program == 'ml':
                metadata['e_slab'] = energies.E_slab_ml
                metadata['e_mol'] = energies.E_mol_ml
            else:
                metadata['e_slab'] = energies.E_slab
                metadata['e_mol'] = energies.E_mol
            return energies
        else:
            # read the energies from the metadata, setting to 0 if were not set before or present from the beginning
            if not metadata.get('e_slab') or not metadata.get('e_mol'):
                logging.warning(f'No slab/molecule energy available in {calc_type} database. ' #pylint: disable=logging-fstring-interpolation
                                'If you calculated slab and mol energies, try to '
                                'refresh the database with "xsorb update -refresh". '
                                'Total energies will be used instead of adsorption energies.')
                metadata['e_slab'] = 0.0
                metadata['e_mol'] = 0.0

            e_slab = metadata['e_slab']
            e_mol = metadata['e_mol']
            if program == 'ml':
                return Energies(E_slab_ml=e_slab, E_mol_ml=e_mol)
            else:
                return Energies(E_slab=e_slab, E_mol=e_mol)


    @staticmethod
    def update_calc_db(calc_type : str,
                            refresh : bool = False,
                            energies : Energies | None = None,
                            mult : float | None = None,
                            store_full_trajectories : bool | None = None,
                            write_csv : bool = True,
                            txt : bool = False,
                            verbose: bool=False) -> None:
        '''
        Update the database with the new results.
        Also update the job status

        Args:
        - calc_type: string with the type of calculation: 'screening'/'relax'/'mlopt', or 'all'
        - refresh: bool to force the update of the database
        - energies: Energies object with slab and mol energies
        - mult: float with the multiplicative factor for the covalent radii to
            determine bonding. Needs to be passed when refreshing the database
            if the value was changed from the settings
        - store_full_trajectories: bool to store the full trajectory in the database
        - write_csv: bool to write the results to a csv file
        - txt: bool to write a txt file instead of a csv file
        - verbose: bool to print messages
        '''

        with ase.db.connect(CALC_DB_NAMES[calc_type]) as db:

            if verbose:
                logging.info('Updating %s database...', calc_type)

            #Get the ids and calc_ids of the (incomplete) calculations to be updated
            selection = 'status!=completed' if not refresh else None
            rows = list(db.select(selection, include_data=True))

            if not rows:
                if verbose:
                    logging.info(
                        'All %s calculations were already completed. ' \
                        'db is already up to date.', calc_type)
                return

            ### read stuff from the database metadata ###
            metadata = db.metadata.copy() # workaround for ase.db bug
            try:
                program = db.metadata['program']
            except KeyError as exc:
                raise RuntimeError(f'No program metadata found in {calc_type} database') from exc


            ### update the metadata with respective values if provided,
            # otherwise read values from metadata
            if mult is not None:
                metadata['mult'] = mult
            else:
                mult = metadata['mult']

            if store_full_trajectories is not None:
                metadata['store_full_trajectories'] = store_full_trajectories
            else:
                store_full_trajectories = metadata['store_full_trajectories']

            #Get the updated slab and mol energies
            updated_energies = Database.get_updated_energies(
                metadata=metadata,
                energies=energies,
                program=program,
                calc_type=calc_type)
            ###############################################


            ### read results and update the calculations ###
            systems = [AdsorptionCalculation.fromdict(row.data.AdsorptionCalculation) \
                       for row in rows]

            #Add the results to the systems
            xsorb.calculations.results.update_calculations_results(
                    systems=systems,
                    program=program,
                    mult=mult,
                    energies=updated_energies,
                    store_full_trajectories=store_full_trajectories,
                    verbose=verbose)
            ################################################

            # Update the database with the new results
            for row, system in zip(rows, systems):
                logging.debug('Updating calculation %s (row %s)', system.calc_info.calc_id, row.id)
                if system.calc_results is not None:
                    db.update(id=row.id,
                            atoms=system.calc_results.atoms,
                            data={'AdsorptionCalculation': system},
                            **system.db_keys())
            db.metadata = metadata # workaround for ase.db bug

        if verbose:
            logging.info('%s database updated.', calc_type)

        if write_csv:
            Database.write_csvfile(txt=txt, verbose=verbose)


    @staticmethod
    def get_structures(calc_ids : list[int] | int | None = None) -> list:
        '''
        Get the rows from the structures database

        Args:
        - calc_ids: list of strings with the ids of the structures to be included,
            or a single integer with the id of the structure to be included.
            If None, all the structures are included

        Returns:
        - list of rows
        '''
        with ase.db.connect(STRUCTURES_DB_NAME) as db:
            rows = list(db.select())
            for row in rows:
                row.__dict__.update({'calc_id': row.id})
        if calc_ids:
            if isinstance(calc_ids, int):
                calc_ids = [calc_ids]
            rows = [row for row in rows if row.id in calc_ids]

        return rows

    #@db_getter
    @staticmethod
    def get_calculations(*,calc_type : str,
                         selection : str | None = None,
                         calc_ids : list[int] | int | None = None,
                         exclude_ids : list[int] | None = None,
                         columns : list[str] | str = 'all',
                         sort_key : str | None = None,
                         include_data : bool = True,
                         update: bool = False,
                         verbose : bool = True) -> list:
        '''
        Get the rows corresponding to the calculations of a given type,
        with the possibility to sort them by a given key

        Args:
        - calc_type: string with the type of calculation
        - selection: string with the selection criteria (e.g. 'status=completed')
        - calc_ids: list of integers with the ids of the calculations to be included,
            or a single integer with the id of the calculation to be included.
            If None, all the calculations are included.
            IMPORTANT: this selects only which ids to include, not the order
        - exclude_ids: list of integers with the ids of the calculations to be excluded
        - columns: list of strings with the columns to be included
        - sort_key: string with the key to sort the rows, e.g. 'energy'
        - include_data: bool to include the data dictionary in the rows
        - update: bool to update the database before getting the calculations
        - verbose: bool to print messages

        Returns:
        - list: list of rows
        '''
        if selection is not None and calc_ids is not None:
            raise ValueError('Cannot use both selection and calc_ids')

        if not Path(CALC_DB_NAMES[calc_type]).exists():
            if verbose:
                logging.warning('Warning: No %s calculations present in the database.', calc_type)
            return []

        #Make sure that the database is up to date
        if update:
            Database.update_calc_db(calc_type, verbose=False)

        with ase.db.connect(CALC_DB_NAMES[calc_type]) as db:
            rows = list(db.select(selection=selection,
                             columns=columns,
                             sort=sort_key,
                             include_data=include_data))
        if calc_ids:
            if isinstance(calc_ids, int):
                calc_ids = [calc_ids]
            rows = [row for row in rows if row.calc_id in calc_ids]
        if exclude_ids:
            rows = [row for row in rows if row.calc_id not in exclude_ids]

        return rows


    @staticmethod
    def remove_calculations(calc_ids: list[int], calc_type: str) -> None:
        '''
        Remove entries from the database for a given list of calc_ids

        Args:
        - calc_ids: list of integers with the ids of the calculations to be removed
        - calc_type: string with the type of calculation (screening, relax, or mlopt)
        '''
        with ase.db.connect(CALC_DB_NAMES[calc_type]) as db:
            for calc_id in calc_ids:
                try:
                    row_id = db.get(f'calc_id={calc_id}', include_data=False).id
                    del db[row_id]
                except KeyError:
                    # If the calculation is not present, do nothing
                    pass


    @staticmethod
    def add_job_id(calc_type : str, calc_id : int, job_id : int) -> None:
        '''
        Add the job id to the corresponding database

        Args:
        - calc_type: string with the type of calculation
        - calc_id: int with the calculation id
        - job_id: int with the job id
        '''
        # we need to do it for each calculation separately,
        # so that if the code crashes or is interrupted,
        # the already submitted jobs are recorded in the database

        with ase.db.connect(CALC_DB_NAMES[calc_type]) as db:
            metadata = db.metadata.copy() # workaround for ase.db bug
            row_id = db.get(f'calc_id={calc_id}', include_data=False).id
            db.update(id=row_id, job_id=job_id)
            db.metadata = metadata # workaround for ase.db bug


    @staticmethod
    def get_all_job_ids() -> list[int]:
        '''
        Get all the job ids from the databases

        Returns:
        - list of integers with the job ids
        '''

        job_ids : list[int] = []
        for calc_type in CALC_DB_NAMES.values():
            if Path(calc_type).exists():
                with ase.db.connect(calc_type) as db:
                    for row in db.select(include_data=False):
                        if row.status!='completed' and row.get('job_id') is not None:
                            job_ids.append(row.job_id)

        return job_ids


    @staticmethod
    def get_adsorption_sites(
            cls_type: AdsorptionSiteCrystal|AdsorptionSiteAmorphous|None=None
            ) -> list[AdsorptionSiteCrystal|AdsorptionSiteAmorphous]:
        '''
        Get the unique adsorption sites from the structures database

        Args:
        - cls_type: class type of the adsorption site to be returned,
            either AdsorptionSiteCrystal or AdsorptionSiteAmorphous

        Returns:
        - list of AdsorptionSiteCrystal or AdsorptionSiteAmorphous objects
        '''

        if not Path(STRUCTURES_DB_NAME).exists():
            logging.debug('Warning: No structures database found. Returning empty list.')
            return []

        with ase.db.connect(STRUCTURES_DB_NAME) as db:
            sites = db.metadata.get('adsorption_sites', [])

        converted_sites = []
        for site in sites:
            if site['__xsorb_objtype__'] == 'AdsorptionSiteCrystal':
                site = AdsorptionSiteCrystal.fromdict(site)
            elif site['__xsorb_objtype__'] == 'AdsorptionSiteAmorphous':
                site = AdsorptionSiteAmorphous.fromdict(site)
            else:
                raise ValueError(f"Unknown site type: {site['__xsorb_objtype__']}")

            if cls_type is None or isinstance(site, cls_type):
                converted_sites.append(site)
            else:
                raise TypeError(f"Expected site of type {cls_type}, got {type(site)}")

        return converted_sites


    #@db_getter
    @staticmethod
    def all_completed(calc_type : str) -> bool:
        '''
        Check if all the calculations in the database are completed

        Args:
        - calc_type: string with the type of calculation

        Returns:
        - bool: True if all the calculations are completed, False otherwise
        '''

        if not Path(CALC_DB_NAMES[calc_type]).exists():
            return False

        #Make sure that the database is up to date
        Database.update_calc_db(calc_type, verbose=False)

        with ase.db.connect(CALC_DB_NAMES[calc_type]) as db:
            return len(list(db.select('status!=completed', include_data=False))) == 0


    @staticmethod
    def write_csvfile(include_results : bool = True,
                      txt : bool = False,
                      sort_txt: bool = False,
                      verbose: bool = False) -> None:
        '''
        Reads the structures database and writes the info to a csv file.
        Is not meant to be called from the CLI, since it is executed
        automatically when adding new structures or updating the calculations.

        It can be called explicitly from the CLI with the command xsorb write_csv
        in case it was accidentally deleted and no calculations are present yet,
        to write the entries of the structures database.
        When results are present, simply call xsorb update to write it again.

        Args:
        - include_results: bool to include the results of the calculations
        - txt: bool to write a txt file instead of a csv file
        - sort_txt: bool to sort the txt output by energy
        - verbose: bool to print messages
        '''

        if verbose: logging.info('Writing results file...') #pylint: disable=multiple-statements

        if not Path(STRUCTURES_DB_NAME).exists():
            raise RuntimeError('Missing structures database. Cannot write csv file.')

        # Get the data from the structures database
        info_dicts : dict[int,dict] = {}
        with ase.db.connect(STRUCTURES_DB_NAME) as db:
            for row in db.select(include_data=False):
                local_dict = {}
                for key in AdsorptionStructure.dataframe_column_names():
                    local_dict.update({key: row.get(key)})
                info_dicts[row.id] = local_dict

        last_calc_e_column_name = None
        if include_results:
            energies_column_names = []

            for calc_type, db_name in CALC_DB_NAMES.items():
                #order is mlopt, screening, relax
                atleast_one_calc = False
                eads_column_name = f'Eads_{calc_type}(eV)'

                if Path(db_name).exists():
                    with ase.db.connect(db_name) as db:
                        for calc_id in info_dicts: #pylint: disable=consider-using-dict-items

                            try:
                                row = db.get(f'calc_id={calc_id}', include_data=False)
                            except KeyError:
                                continue #calc_id not present in this db

                            eads = row.get('adsorption_energy') #can be a float or None
                            if eads is not None:
                                atleast_one_calc = True
                                eads = f'{eads:.3f}'
                                if row.get('status') == 'scf_nonconverged':
                                    eads += '**'
                                    if verbose:
                                        logging.warning('Warning! %s %s '
                                          'failed to reach SCF convergence in the last step. '
                                            'The energy will be marked with **', calc_type, calc_id)
                                elif row.get('status') == 'incomplete':
                                    eads += '*'
                                    if verbose:
                                        logging.warning('Warning! %s %s '
                                          'has not reached final configuration. '
                                            'The energy will be marked with a *', calc_type, calc_id)

                            info_dicts[calc_id].update({eads_column_name: eads})

                            info_dicts[calc_id].update({f'bonds_{calc_type}': row.get('bonds')})

                            info_dicts[calc_id].update({f'final_dz_{calc_type}': row.get("final_dz")})


                    if atleast_one_calc:
                        last_calc_e_column_name = eads_column_name
                        energies_column_names.append(eads_column_name)

        # Write csv file
        df = pd.DataFrame.from_dict(info_dicts, orient='index')

        # change index to calc_id
        df.index.name = 'calc_id'
        df.reset_index(inplace=True)

        if last_calc_e_column_name and include_results and txt and sort_txt: #sort by energy column
            print(f'Sorting results by {last_calc_e_column_name}')
            df.sort_values(by=last_calc_e_column_name, inplace=True, ascending=False)

        if txt:
            with open('results.txt', 'w') as f:
                f.write(df.to_string(index=False, float_format='%.2f', na_rep='-', col_space=2))
        else:
            df.to_csv('results.csv', index=False, float_format='%.3f')

        if verbose: logging.info('Results file written.') #pylint: disable=multiple-statements


def manual_update_calculations(calc_type : str,
                               refresh : bool = False,
                               txt : bool = False,
                               sort_txt: bool = False) -> None:
    '''
    Manually update the calculations in the database and write the results file.
    This function is meant to be called from the CLI

    Args:
    - calc_type: string with the type of calculation to update
    - refresh: bool to force the update of the database, re-reading all the output files
      and updating the bonding status.
    - txt: bool to write a txt file instead of a csv file
    - sort_txt: bool to sort the txt output by energy
    '''

    if refresh:
        logging.info('Re-reading the output files, updating slab and mol energies, '
                'the radii mult factor, and recalculating the bonding status...')
        from xsorb.settings import Settings # pylint: disable=import-outside-toplevel
        settings = Settings(read_energies_dft=(calc_type != 'mlopt'),
                            read_energies_ml=(calc_type == 'mlopt' or calc_type == 'all'))
        mult=settings.structure.molecule.radius_scale_factor
        energies = settings.energies
        store_full_trajectories = settings.database.store_full_trajectories
    else:
        mult = None
        energies = None
        store_full_trajectories = None


    dbs_to_update = [calc_type] if calc_type != 'all' else CALC_DB_NAMES.keys()
    for calc_type in dbs_to_update:
        if Path(CALC_DB_NAMES[calc_type]).exists():
            Database.update_calc_db(
                calc_type=calc_type,
                refresh=refresh,
                energies=energies,
                mult=mult,
                store_full_trajectories=store_full_trajectories,
                write_csv=False,
                verbose=True
            )

    # Write the results file only once at the end
    Database.write_csvfile(txt=txt, sort_txt=sort_txt, verbose=True)
