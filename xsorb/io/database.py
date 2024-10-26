#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Enrico Pedretti

'''
Module to handle the ase database

Four databases are used:
- structures.db: Database with the initial structures
- screening.db: Contains the results of the screening calculations
- relaxations.db: Contains the results of the relaxations
- ml_opt.db: Contains the results of the machine learning optimization


The calculation databases have the following columns:
- id (int): unique id of the row
- calc_id (int): unique id of the calculation
- atoms (Atoms): Atoms object with the structure
- constraints (list): list of constraints applied to the structure. May not be there
- energy (float): total energy of the structure. May not be there
- status (str): status of the calculation (incomplete, completed)
- scf_nonconverged (bool): True if the SCF did not converge
- job_id (int): id of the job in the scheduler
- job_status (str): status of the job in the scheduler
- in_file_path (str): path to the input file
- out_file_path (str): path to the output file
- log_file_path (str): path to the log file
- the "site_info", "x", "y", "z", "initial_dz", "xrot", "yrot", "zrot" keys
    from AdsorptionStructure.dataframe_column_names
- adsorption_energy (float): adsorption energy of the structure
- bonds (list): list of bonds in the structure
- final_dz (float): final vertical distance between the reference atom of the molecule
    and the adsorption site
- data (dict): dictionary with additional data, namely
    -the AdsorptionStructure object,
    -the 'trajectory'
    -the 'adsorption_energy_evolution'
    and the AdsorptionStructure.additional_data_arrays():
    - slab_indices
    - mol_indices
    - custom_labels

Each database also has the following metadata:
- program (str): name of the program used for the calculations
- mult (float): multiplicative factor for the covalent radii to determine bonding
    The mult factor is updated from settings only when refreshing the database
- total_e_slab_mol (float): total energy of the isolated molecule and slab. May be None

'''
from __future__ import annotations
from typing import TYPE_CHECKING
from pathlib import Path

import pandas as pd
import ase.db
import ase.db.core

import xsorb.calculations.results
if TYPE_CHECKING:
    from xsorb.io.inputs import WrittenSystem
    from xsorb.structures.generation import AdsorptionStructure

DATAFRAME_COLUMNS_NAMES = ("site", "site_info", "mol_atom", "initial_dz", "xrot", "yrot", "zrot")


#TODO: check that the atoms object stored in the database also contains custom_labels. Otherwise,
# make sure they are returned when necessary from the data key

class Database:
    '''
    Collection of static functions to read/write generated structures and calculations
    results to ase db
    '''

    calc_types = {
        'ml_opt': 'ml_opt.db',
        'screening': 'screening.db',
        'relax': 'relaxations.db',
    }


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
        - list of integers with the calc_ids of the structures
        '''

        # Write the adsorption structures to the database
        # excluding those that are already present
        calc_ids = []
        with ase.db.connect('structures.db') as db:
            for ads_struct in adsorption_structures:
                already_present = False
                for row in db.select(include_data=False):
                    if ads_struct.atoms == row.toatoms():
                        already_present = True
                        calc_ids.append(row.id)
                        break
                if not already_present:
                    calc_id = db.write(ads_struct.atoms,
                                       kwargs=ads_struct.to_info_dict(),
                                       data=ads_struct.additional_data_arrays())
                    calc_ids.append(calc_id)

        if write_csv:
            Database.write_csvfile(include_results=False)

        return calc_ids


    @staticmethod
    def add_calculations(systems : list [WrittenSystem],
                         program : str,
                         mult : float,
                         total_e_slab_mol : float | None,
                         calc_type : str) -> None:
        '''
        Write new calculations to corresponding database

        Args:
        - systems: list of WrittenSystem objects, each containing
            'calc_id', 'adsorption_structure', 'in_file_path', 'out_file_path'
        - program: string with the name of the program used for the calculations:
            'vasp', 'espresso', 'ml'
        - mult: float with the multiplicative factor for the covalent radii to
            determine bonding
        - total_e_slab_mol: float with the total energy of the isolated molecule and slab
        - calc_type: screening, relax or ml_opt
        '''

        # Write the adsorption structures to the corresponding database
        with ase.db.connect(Database.calc_types[calc_type]) as db:

            db.metadata = {'program': program,
                           'mult': mult,
                           'total_e_slab_mol': total_e_slab_mol}

            for system in systems:
                try:
                    #if the calculation is already present in the database, remove it
                    del db[db.get(f'calc_id={system.calc_id}').id]
                except KeyError:
                    #if the calculation is not present, do nothing
                    pass

                #write the new calculation in any case
                ads_struct : AdsorptionStructure = system.adsorption_structure
                data=ads_struct.additional_data_arrays()
                data.update({'adsorption_structure': ads_struct})
                db.write(ads_struct.atoms,
                        calc_id=system.calc_id,
                        kwargs=ads_struct.to_info_dict(),
                        status='incomplete',
                        in_file_path=system.in_file_path,
                        out_file_path=system.out_file_path,
                        log_file_path=system.log_file_path,
                        data=data)


    @staticmethod
    def update_calculations(calc_type : str,
                            refresh : bool = False,
                            mult : float | None = None,
                            write_csv : bool = True,
                            txt : bool = False) -> None:
        '''
        Update the database with the new results.
        Also update the job status

        Args:
        - calc_type: string with the type of calculation: 'screening'/'relax'/'ml_opt', or 'all'
        - refresh: bool to force the update of the database
        - mult: float with the multiplicative factor for the covalent radii to
            determine bonding. Needs to be passed when refreshing the database
            if the value was changed from the settings
        - write_csv: bool to write the results to a csv file
        - txt: bool to write a txt file instead of a csv file
        '''

        if calc_type == 'all':
            for ctype, db_name in Database.calc_types.items():
                if Path(db_name).exists():
                    Database.update_calculations(ctype, refresh, mult, write_csv=False)

            if write_csv: #only write the csv file once
                Database.write_csvfile(txt=txt)
            return

        with ase.db.connect(Database.calc_types[calc_type]) as db:
            #Get the ids and calc_ids of the (incomplete) calculations to be updated
            selection = 'status=incomplete' if not refresh else None
            rows = db.select(selection, include_data=True)

            row_ids = [row.id for row in rows]

            program = db.metadata.get('program')
            if mult is not None:
                db.metadata['mult'] = mult
            else:
                mult = db.metadata.get('mult')
            total_e_slab_mol = db.metadata.get('total_e_slab_mol')
            systems = [WrittenSystem(calc_id='',
                                     adsorption_structure=None,
                                     in_file_path=row.get('in_file_path'),
                                     out_file_path=row.get('out_file_path'),
                                     log_file_path=row.get('log_file_path'),
                                     job_id=row.get('job_id'),
                                     adsite_z=row.data.adsorption_structure.adsite.coords[2],
                                     mol_ref_idx=row.data.adsorption_structure.mol_rot.mol_atom)
                                        for row in rows]

            #Get the results of the calculations
            results = xsorb.calculations.results.get_calculations_results(
                    systems=systems,
                    program=program,
                    mol_indices=rows[0].data.mol_indices,
                    mult=mult,
                    total_e_slab_mol=total_e_slab_mol)

            for row_id, result in zip(row_ids, results):
                if results is not None:
                    db.update(id=row_id,
                              atoms=result.atoms,
                              status=result.status,
                              scf_nonconverged=result.scf_nonconverged,
                              adsorption_energy=result.adsorption_energy,
                              bonds=result.bonds,
                              final_dz=result.final_dz,
                              job_status=result.job_status,
                              data={'trajectory': result.trajectory,
                                    'adsorption_energy_evolution': result.adsorption_energy_evol})
                else:
                    pass
                    #this is probably already printed somewhere else. Check
                    #print(f'Warning: Calculation {row.calc_id} cannot be updated.')

        if write_csv:
            Database.write_csvfile(txt=txt)

    @staticmethod
    def get_structures(calc_ids : list[int] | int | None = None) -> list:
        '''
        Get the structures from the structures database

        Args:
        - calc_ids: list of integers with the ids of the structures to be included,
            or a single integer with the id of the structure to be included.
            If None, all the structures are included

        Returns:
        - list of rows
        '''
        with ase.db.connect('structures.db') as db:
            rows = db.select(include_data=False)
        if calc_ids:
            if isinstance(calc_ids, int):
                calc_ids = [calc_ids]
            rows = [row for row in rows if row.id in calc_ids]

        return rows

    #@db_getter
    @staticmethod
    def get_calculations(calc_type : str,
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

        if not Path(Database.calc_types[calc_type]).exists():
            print(f'Warning: No {calc_type} calculations present in the database.')
            return []

        #Make sure that the database is up to date
        Database.update_calculations(calc_type)

        with ase.db.connect(Database.calc_types[calc_type]) as db:
            rows = db.select(selection=selection,
                             columns=columns,
                             sort=sort_key,
                             include_data=include_data)
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
        - calc_type: string with the type of calculation (screening, relax, or ml_opt)
        '''
        with ase.db.connect(Database.calc_types[calc_type]) as db:
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
        with ase.db.connect(Database.calc_types[calc_type]) as db:
            for calc_id, job_id in zip(calc_ids, job_ids):
                row_id = db.get(f'calc_id={calc_id}', include_data=False).id
                db.update(id=row_id, job_id=job_id, job_status='submitted')


    @staticmethod
    def get_all_job_ids() -> list:
        '''
        Get all the job ids from the databases

        Returns:
        - list of integers with the job ids
        '''
        job_ids = []
        calc_types = Database.calc_types.copy()

        for calc_type in calc_types.values():
            if Path(calc_type).exists():
                with ase.db.connect(calc_type) as db:
                    for row in db.select(include_data=False):
                        if row.status=='incomplete' and row.get('job_id') is not None:
                            job_ids.append(row.job_id)

        return job_ids


    # @staticmethod
    # def get_adsorption_sites(calc_type : str) -> list:
    #     '''
    #     Get the unique adsorption sites from the database

    #     Args:
    #     - calc_type: string with the type of calculation

    #     Returns:
    #     - list of strings with the adsorption sites
    #     '''
    #     with ase.db.connect(Database.calc_types[calc_type]) as db:
    #         return list(set(row.get('site') for row in db.select(include_data=False)))

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

        if not Path(Database.calc_types[calc_type]).exists():
            return False

        #Make sure that the database is up to date
        Database.update_calculations(calc_type)

        with ase.db.connect(Database.calc_types[calc_type]) as db:
            return len(db.select('status=incomplete', include_data=False)) == 0


    @staticmethod
    def write_csvfile(include_results : bool = True, txt : bool = False) -> None:
        '''
        Reads the structures.db database and writes the info to a csv file.
        Is not meant to be called from the CLI, since it is executed
        automatically when adding new structures or updating the calculations.

        It can be called explicitly from the CLI with the command xsorb write_csv
        in case it was accidentally deleted and no calculations are present yet,
        to write the entries of the structures.db.
        When results are present, simply call xsorb update to write it again.

        Args:
        - include_results: bool to include the results of the calculations
        - txt: bool to write a txt file instead of a csv file
        '''

        print('Writing results file...')

        if not Path('structures.db').exists():
            raise RuntimeError('Missing structures.db database. Cannot write csv file.')

        # Get the data from the database
        info_dicts : list[dict] = []
        with ase.db.connect('structures.db') as db:
            for i, row in db.select(include_data=False):
                info_dict = {'calc_id': row.id}
                for key in DATAFRAME_COLUMNS_NAMES:
                    info_dict.update({key: row.get(key)})
                info_dicts.append(info_dict)

        last_calc_e_column_name = None
        if include_results:
            for calc_type, db_name in Database.calc_types.items():
                #order is ml_opt, screening, relax
                if Path(db_name).exists():
                    with ase.db.connect(db_name) as db:
                        for i, info_dict in enumerate(info_dicts):
                            try:
                                row = db.get(f'calc_id={info_dict["calc_id"]}', include_data=False)
                            except KeyError:
                                continue

                            eads = row.get('adsorption_energy') #can be a float or None
                            if eads is not None:
                                eads = f'{eads:.3f}'
                                if row.get('scf_nonconverged'):
                                    eads += '**'
                                    print(f'Warning! {calc_type} {info_dict["calc_id"]} '\
                                          'failed to reach SCF convergence in the last step. '\
                                            'The energy will be marked with **')
                                elif row.get('status') != 'completed':
                                    eads += '*'
                                    print(f'Warning! {calc_type} {info_dict["calc_id"]} '\
                                          'has not reached final configuration. '\
                                            'The energy will be marked with a *')
                            info_dicts[i].update({f'Eads_{calc_type[:3]}(eV)': eads})
                            info_dicts[i].update({'bonds': row.get('bonds')})
                            info_dicts[i].update({'final_dz': row.get('final_dz')})

                    last_calc_e_column_name = f'Eads_{calc_type[:3]}(eV)'

        # Write csv file
        df = pd.DataFrame(columns=['calc_id'] + DATAFRAME_COLUMNS_NAMES)
        for info_dict in info_dicts:
            df_line = pd.Series(info_dict)
            df.loc[i] = df_line

        if include_results and txt: #sort by energy column
            df.sort_values(by=last_calc_e_column_name)

        if txt:
            df.to_csv('results.txt', sep='\t', index=False)
        else:
            df.to_csv('results.csv', index=False)

        print('Results file written.')
