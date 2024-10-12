'''
Module for to handle the ase database

Four databases are used:
- structures.db: Database with the initial structures
- screening.db: Contains the results of the screening calculations
- relaxations.db: Contains the results of the relaxations
- ml_opt.db: Contains the results of the machine learning optimization

'''

import ase.db
from ase import Atoms

from xsorb.structures import AdsorptionStructure
from xsorb.dft_codes.definitions import IN_FILE_PATHS, OUT_FILE_PATHS, LOG_FILE_PATHS

#TODO: add option to remove entries from the database.
#this is useful if the user is unstisfied with some of the generated structures.



#DESIGN CHOICE: the database has ABSOLUTE PRIORITY over files.
#When a new structure is created, it is ALWAYS written to the structures.db
#If the user wants to play around a bit at the beginning to find the optimal parameters, they
#can remove the structures.db file and start from scratch.
#The other databases are updated only when not in generation mode.
#The input files will be overwritten unless they are in their respective database.
#In that case, the user will be asked if they want to overwrite the file or not.
#If entry is in database and user wants to overwrite, check if the folder is present
#and remove it before writing the new files.


class Database:
    '''
    Database class to read/write generated structures and calculations
    results to ase db
    '''

    calc_types = {
        'screening': 'screening.db',
        'ml_opt': 'ml_opt.db',
        'relax': 'relaxations.db',
    }


    @staticmethod
    def add_structures(adsorption_structures : list [AdsorptionStructure]) -> None:
        '''
        Adds the generated adsorption structures to the structures database,
        returning the corresponding calc_ids

        Args:
        - adsorption_structures: list of AdsorptionStructure objects

        Returns:
        - list of integers with the calc_ids of the structures
        '''

        # Write the adsorption structures to the database
        # excluding those that are already present
        calc_ids = []
        with ase.db.connect('structures.db') as db:
            for ads_struct in adsorption_structures:
                already_present = False
                for row in db.select():
                    if ads_struct.atoms == row.toatoms():
                        already_present = True
                        calc_ids.append(row.id)
                        break
                if not already_present:
                    calc_id = db.write(ads_struct.atoms,
                                       kwargs=ads_struct.to_info_dict(),
                                       data=ads_struct.additional_data_arrays(),)
                    calc_ids.append(calc_id)

        return calc_ids




    @staticmethod
    def add_calculations(systems : list [dict],
                         calc_type : str = None) -> None:
        '''
        Write new calculations to corresponding database

        Args:
        - systems: list of dictionaries, each containing
            'calc_id', 'atoms', 'in_file_path', 'out_file_path'
        - calc_type: screening, relax or ml_opt
        '''

        # Write the adsorption structures to the corresponding database
        with ase.db.connect(Database.calc_types[calc_type]) as db:
            for system in systems:
                try:
                    db.get_atoms(f'calc_id={system.get("calc_id")}')
                except KeyError:
                    already_present = False
                if not already_present:
                    db.write(ads_struct.atoms,
                            calc_id=calc_id,
                            kwargs=ads_struct.to_info_dict(),
                            in_filename=IN_FILE_PATHS[calc_type][program].format(calc_id),
                            out_filename=OUT_FILE_PATHS[calc_type][program].format(calc_id),
                            data=ads_struct.additional_data_arrays(),
                               )


    @staticmethod
    def update(atoms_list : list[Atoms],
               calc_ids : list[int],
               convergence_info : dict,
               calc_type : str) -> None:
        '''
        Update the database with the new results

        Args:
        - atoms_list: list of Atoms objects with the new results
        - calc_ids: list of integers with the ids of the calculations
        - convergence_info: dictionary with the convergence information
        - calc_type: string with the type of calculation
        '''
        with ase.db.connect(Database.calc_types[calc_type]) as db:
            for atoms, calc_id in zip(atoms_list, calc_ids):
                db.update(calc_id, atoms=atoms, **convergence_info)


    @staticmethod
    def get(calc_type : str, id :int = None) -> list[Atoms]:
        '''
        Get the atoms objects from the database

        Args:
        - calc_type: string with the type of calculation
        - kwargs: keyword arguments to filter the results

        Returns:
        - list of Atoms objects
        '''
        db = ase.db.connect(Database.calc_types[calc_type])


