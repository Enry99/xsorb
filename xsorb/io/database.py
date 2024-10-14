'''
Module to handle the ase database

Four databases are used:
- structures.db: Database with the initial structures
- screening.db: Contains the results of the screening calculations
- relaxations.db: Contains the results of the relaxations
- ml_opt.db: Contains the results of the machine learning optimization

'''

import numpy as np
import pandas as pd
import ase.db
import ase.db.core

from xsorb.structures import AdsorptionStructure
from xsorb.io.jobs import get_running_jobs


#TODO: convergence_info must contain a {'status': 'completed'/'incomplete'} key

class Database:
    '''
    Collection of functions to read/write generated structures and calculations
    results to ase db
    '''

    calc_types = {
        'screening': 'screening.db',
        'ml_opt': 'ml_opt.db',
        'relax': 'relaxations.db',
    }


    @staticmethod
    def add_structures(adsorption_structures : list [AdsorptionStructure],
                       write_csv : bool = True) -> None:
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
            Database.write_csvfile()

        return calc_ids


    @staticmethod
    def add_calculations(systems : list [dict],
                         calc_type : str = None) -> None:
        '''
        Write new calculations to corresponding database

        Args:
        - systems: list of dictionaries, each containing
            'calc_id', 'adsorption_structure', 'in_file_path', 'out_file_path'
        - calc_type: screening, relax or ml_opt
        '''

        # Write the adsorption structures to the corresponding database
        with ase.db.connect(Database.calc_types[calc_type]) as db:
            for system in systems:
                try:
                    #if the calculation is already present in the database, remove it
                    del db[db.get(f'calc_id={system.get("calc_id")}').id]
                except KeyError:
                    #if the calculation is not present, do nothing
                    pass

                #write the new calculation in any case
                ads_struct : AdsorptionStructure = system['adsorption_structure']
                data=ads_struct.additional_data_arrays()
                data.update({'adsorption_structure': ads_struct})
                db.write(ads_struct.atoms,
                        calc_id=system['calc_id'],
                        kwargs=ads_struct.to_info_dict(),
                        status='incomplete',
                        in_filename=system['in_file_path'],
                        out_filename=system['out_file_path'],
                        data=data)


    @staticmethod
    def update_calculations(calc_type : str):
        '''
        Update the database with the new results.
        Also update the job status

        Args:
        - calc_type: string with the type of calculation

        '''
        #TODO: write here the code to retrieve the data from the output files

        #- atoms_list: list of trajectories
        #- calc_ids: list of integers with the ids of the calculations
        #- convergence_info: list of dictionaries with the convergence information

        # with ase.db.connect(Database.calc_types[calc_type]) as db:
        #     for traj, calc_id, convergence_info in zip(traj_list, calc_ids, convergence_info_list):
        #         row_id = db.get(f'calc_id={calc_id}').id
        #         db.update(id=row_id, atoms=traj[-1], data={'trajectory': traj}, **convergence_info)

        pass

    #@db_getter
    @staticmethod
    def get_calculations(calc_type : str,
                         selection : str | None = None,
                         calc_ids : list[int] | None = None,
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
        - calc_ids: list of integers with the ids of the calculations to be included
        - exclude_ids: list of integers with the ids of the calculations to be excluded
        - columns: list of strings with the columns to be included
        - sort_key: string with the key to sort the rows, e.g. 'energy'

        Returns:
        - list: list of rows
        '''
        if selection is not None and calc_ids is not None:
            raise ValueError('Cannot use both selection and calc_ids')

        #Make sure that the database is up to date
        Database.update_calculations(calc_type)

        with ase.db.connect(Database.calc_types[calc_type]) as db:
            rows = db.select(selection=selection,
                             columns=columns,
                             sort=sort_key,
                             include_data=include_data)
        if calc_ids:
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
    def write_csvfile():
        '''
        Reads the structures.db database and writes the info to a csv file
        '''

        df = pd.DataFrame(columns=['calc_id'] + AdsorptionStructure.dataframe_column_names)

        with ase.db.connect('structures.db') as db:
            for i, row in db.select(include_data=False):
                info_dict = {'calc_id': row.id}
                for key in AdsorptionStructure.dataframe_column_names:
                    info_dict.update(key, row.get(key))
                df_line = pd.Series(info_dict)
                df.loc[i] = df_line

        df.to_csv('results.csv', index=False)


    @staticmethod
    def add_job_ids(calc_type : str, calc_ids : list[int], job_ids : list[int]) -> None:
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
    def update_job_status(calc_type : str) -> None:
        '''
        Go over the database and update the job status,
        for the calculations that are not already marked as completed,
        by checking if the job_id is present in the scheduler output
        '''

        running_jobs = get_running_jobs()

        with ase.db.connect(Database.calc_types[calc_type]) as db:
            for row in db.select(include_data=False):
                if row.status != 'completed':
                    job_id = row.job_id
                    if job_id in running_jobs:
                        job_status='running'
                    else:
                        job_status='terminated'
                    db.update(id=row.id, job_status=job_status)


    @staticmethod
    def get_all_job_ids() -> list:
        '''
        Get all the job ids from the databases

        Returns:
        - list of integers with the job ids
        '''
        job_ids = []
        calc_types = Database.calc_types.copy()
        calc_types.pop('structures')
        for calc_type in calc_types.values():
            with ase.db.connect(calc_type) as db:
                for row in db.select(include_data=False):
                    if row.status=='incomplete' and row.get('job_id') is not None:
                        job_ids.append(row.job_id)

        return job_ids


    @staticmethod
    def get_adsorption_sites(calc_type : str) -> list:
        '''
        Get the unique adsorption sites from the database

        Args:
        - calc_type: string with the type of calculation

        Returns:
        - list of strings with the adsorption sites
        '''
        with ase.db.connect(Database.calc_types[calc_type]) as db:
            return list(set(row.get('site') for row in db.select(include_data=False)))

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

        #Make sure that the database is up to date
        Database.update_calculations(calc_type)

        with ase.db.connect(Database.calc_types[calc_type]) as db:
            return len(db.select('status=incomplete', include_data=False)) == 0
