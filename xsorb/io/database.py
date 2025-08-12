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
- "site", "site_info", "xrot", "yrot", "zrot", "mol_atom", "coords", "initial_dz"

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
- total_e_slab_mol (float): total energy of the isolated molecule and slab. May be None
while the structure database contains:
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
from xsorb.adsorptiondata.adsorptionstructure import AdsorptionSiteCrystal, AdsorptionSiteAmorphous
from xsorb.adsorptiondata import AdsorptionCalculation
from xsorb.adsorptiondata import AdsorptionStructure


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
                if ads_struct.adsite not in adsorption_sites:
                    adsorption_sites.append(ads_struct.adsite)

            # Update the metadata with the adsorption sites
            db.metadata['adsorption_sites'] = [site for site in adsorption_sites]

        if write_csv:
            Database.write_csvfile(include_results=False)

        return calc_ids


    @staticmethod
    def add_calculations(systems : list [AdsorptionCalculation],
                         program : str,
                         mult : float,
                         total_e_slab_mol : float | None,
                         calc_type : str) -> None:
        '''
        Write new calculations to corresponding database

        Args:
        - systems: list of CalculationInfo objects, each containing
            'calc_id', 'adsorption_structure', 'in_file_path', 'out_file_path'
        - program: string with the name of the program used for the calculations:
            'vasp', 'espresso', 'ml'
        - mult: float with the multiplicative factor for the covalent radii to
            determine bonding
        - total_e_slab_mol: float with the total energy of the isolated molecule and slab
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
                        **ads_struct.db_keys())

            db.metadata = {'program': program,
                           'mult': mult,
                           'total_e_slab_mol': total_e_slab_mol}


    @staticmethod
    def update_calc_db(calc_type : str,
                            refresh : bool = False,
                            total_e_slab_mol_dft : float | None = None,
                            total_e_slab_mol_ml : float | None = None,
                            mult : float | None = None,
                            write_csv : bool = True,
                            txt : bool = False,
                            verbose: bool=False) -> None:
        '''
        Update the database with the new results.
        Also update the job status

        Args:
        - calc_type: string with the type of calculation: 'screening'/'relax'/'mlopt', or 'all'
        - refresh: bool to force the update of the database
        - total_e_slab_mol_dft: float with the dft total energy of the isolated molecule and slab.
        - total_e_slab_mol_ml: float with total energy of the isolated molecule and slab (for ML)
        - mult: float with the multiplicative factor for the covalent radii to
            determine bonding. Needs to be passed when refreshing the database
            if the value was changed from the settings
        - write_csv: bool to write the results to a csv file
        - txt: bool to write a txt file instead of a csv file
        - verbose: bool to print messages
        '''

        with ase.db.connect(CALC_DB_NAMES[calc_type]) as db:
            #Get the ids and calc_ids of the (incomplete) calculations to be updated
            selection = 'status!=completed' if not refresh else None
            rows = list(db.select(selection, include_data=True))

            ### read stuff from the database metadata ###
            try:
                program = db.metadata.get('program')
            except KeyError as exc:
                raise RuntimeError(f'No program metadata found in {calc_type} database') from exc


            if mult is not None:
                db.metadata['mult'] = mult
            else:
                try:
                    mult = db.metadata.get('mult')
                except KeyError as exc:
                    raise RuntimeError(f'No mult metadata found in {calc_type} database') from exc

            if calc_type == 'mlopt':
                if total_e_slab_mol_ml is None: # read from metadata
                    try:
                        total_e_slab_mol = db.metadata.get('total_e_slab_mol')
                    except KeyError as exc:
                        raise RuntimeError('No total_e_slab_mol metadata found in '\
                                           f'{calc_type} database') from exc
                else: # use the provided value, and update the metadata
                    total_e_slab_mol = total_e_slab_mol_ml
                    db.metadata['total_e_slab_mol'] = total_e_slab_mol
            else:
                if total_e_slab_mol_dft is None: # read from metadata
                    try:
                        total_e_slab_mol = db.metadata.get('total_e_slab_mol')
                    except KeyError as exc:
                        raise RuntimeError('No total_e_slab_mol metadata found in '\
                                           f'{calc_type} database') from exc
                else: # use the provided value, and update the metadata
                    total_e_slab_mol = total_e_slab_mol_dft
                    db.metadata['total_e_slab_mol'] = total_e_slab_mol
            ### end of reading metadata ###


            systems = [AdsorptionCalculation.fromdict(row.data.AdsorptionCalculation) \
                       for row in rows]

            #Add the results to the systems
            xsorb.calculations.results.update_calculations_results(
                    systems=systems,
                    program=program,
                    mult=mult,
                    total_e_slab_mol=total_e_slab_mol,
                    verbose=verbose)

        if verbose:
            logging.info(f'{calc_type} database updated.')

        if write_csv:
            Database.write_csvfile(txt=txt, verbose=verbose)


    @staticmethod
    def update_all_calc_dbs(*,refresh : bool = False,
                            total_e_slab_mol_dft : float | None = None,
                            total_e_slab_mol_ml : float | None = None,
                            mult : float | None = None,
                            write_csv : bool = True,
                            txt : bool = False,
                            verbose: bool=False) -> None:
        '''
        Update all the calculation databases

        Args:
        - refresh: bool to force the update of the database
        - total_e_slab_mol_dft: float with the dft total energy of the isolated molecule and slab.
        - total_e_slab_mol_ml: float with total energy of the isolated molecule and slab (for ML)
        - mult: float with the multiplicative factor for the covalent radii to
            determine bonding. Needs to be passed when refreshing the database
            if the value was changed from the settings
        - write_csv: bool to write the results to a csv file
        - txt: bool to write a txt file instead of a csv file
        - verbose: bool to print messages
        '''


        if verbose and refresh:
            logging.info('Re-reading the output files, updating e_slab_mol, '\
                   'the radii mult factor, and recalculating the bonding status...')


        for calc_type in CALC_DB_NAMES:
            Database.update_calc_db(calc_type=calc_type,
                                    refresh=refresh,
                                    total_e_slab_mol_dft=total_e_slab_mol_dft,
                                    total_e_slab_mol_ml=total_e_slab_mol_ml,
                                    mult=mult,
                                    write_csv=False,
                                    txt=txt,
                                    verbose=verbose)

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
                         include_data : bool = True) -> list:
        '''
        Get the rows corresponding to the calculations of a given type,
        with the possibility to sort them by a given key

        Args:
        - calc_type: string with the type of calculation
        - selection: string with the selection criteria (e.g. 'status=completed')
        - calc_ids: list of integers with the ids of the calculations to be included,
            or a single integer with the id of the calculation to be included.
            If None, all the calculations are included
        - exclude_ids: list of integers with the ids of the calculations to be excluded
        - columns: list of strings with the columns to be included
        - sort_key: string with the key to sort the rows, e.g. 'energy'

        Returns:
        - list: list of rows
        '''
        if selection is not None and calc_ids is not None:
            raise ValueError('Cannot use both selection and calc_ids')

        if not Path(CALC_DB_NAMES[calc_type]).exists():
            logging.warning(f'Warning: No {calc_type} calculations present in the database.')
            return []

        #Make sure that the database is up to date
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
    def add_job_ids(calc_type : str, calc_ids : list[int], job_ids : list[str]) -> None:
        '''
        Add the job ids to the corresponding database

        Args:
        - calc_type: string with the type of calculation
        - calc_ids: list of integers with the calculation ids
        - job_ids: list of integers with the job ids
        '''
        with ase.db.connect(CALC_DB_NAMES[calc_type]) as db:
            for calc_id, job_id in zip(calc_ids, job_ids):
                row_id = db.get(f'calc_id={calc_id}', include_data=False).id
                db.update(id=row_id, job_id=job_id)


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
    def get_adsorption_sites(cls_type=None) -> list[AdsorptionSiteCrystal|AdsorptionSiteAmorphous]:
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
                eads_column_name = f'Eads_{calc_type[:3]}(eV)'

                if Path(db_name).exists():
                    with ase.db.connect(db_name) as db:
                        for calc_id in info_dicts: #pylint: disable=consider-using-dict-items
                            try:
                                row = db.get(f'calc_id={calc_id}', include_data=False)
                            except KeyError:
                                continue


                            eads = row.get('adsorption_energy') #can be a float or None
                            if eads is not None:
                                eads = f'{eads:.3f}'
                                if row.get('status') == 'scf_nonconverged':
                                    eads += '**'
                                    if verbose:
                                        logging.warning(f'Warning! {calc_type} {calc_id} '\
                                          'failed to reach SCF convergence in the last step. '\
                                            'The energy will be marked with **')
                                elif row.get('status') == 'incomplete':
                                    eads += '*'
                                    if verbose:
                                        logging.warning(f'Warning! {calc_type} {calc_id} '\
                                          'has not reached final configuration. '\
                                            'The energy will be marked with a *')
                            info_dicts[calc_id].update({eads_column_name: eads})

                            if row.get('bonds'):
                                info_dicts[calc_id].update({'bonds': row.get('bonds')})

                            if row.get('final_dz'):
                                info_dicts[calc_id].update({'final_dz': row.get('final_dz')})

                            atleast_one_calc = True

                    if atleast_one_calc:
                        last_calc_e_column_name = eads_column_name
                        energies_column_names.append(eads_column_name)

        # Write csv file
        df = pd.DataFrame.from_dict(info_dicts, orient='index')

        if last_calc_e_column_name and include_results and txt: #sort by energy column
            df.sort_values(by=last_calc_e_column_name)

        if txt:
            df.to_csv('results.txt', sep='\t', index=False)
        else:
            df.to_csv('results.csv', index=False)

        if verbose: logging.info('Results file written.') #pylint: disable=multiple-statements


def manual_update_calculations(calc_type : str,
                               refresh : bool = False,
                               txt : bool = False) -> None:
    '''
    Manually update the calculations in the database.
    This function is meant to be called from the CLI, to update the database

    Args:
    - calc_type: string with the type of calculation to update
    - refresh: bool to force the update of the database, re-reading all the output files
      and updating the bonding status.
    - txt: bool to write a txt file instead of a csv file
    '''

    if refresh:
        from xsorb.settings import Settings # pylint: disable=import-outside-toplevel
        settings = Settings(read_energies=True)
        total_e_slab_mol = settings.total_e_slab_mol
        total_e_slab_mol_ml = settings.total_e_slab_mol_ml
        mult=settings.structure.molecule.radius_scale_factor
    else:
        mult = None
        total_e_slab_mol = None
        total_e_slab_mol_ml = None


    Database.update_calc_db(
        calc_type=calc_type,
        refresh=refresh,
        total_e_slab_mol_dft=total_e_slab_mol,
        total_e_slab_mol_ml=total_e_slab_mol_ml,
        mult=mult,
        txt=txt,
        verbose=True
    )
