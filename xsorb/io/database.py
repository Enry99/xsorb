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

class Database:
    '''
    Database class to read/write generated structures and calculations
    results to ase db
    '''

    calc_types = {
        'structures': 'structures.db',
        'screening': 'screening.db',
        'relaxation': 'relaxations.db',
        'ml_opt': 'ml_opt.db'
    }

    @staticmethod
    def write_structures(adsorption_structures : list [AdsorptionStructure],
                         filenames : list[str],
                         calc_type : str = None) -> None:
        '''
        Write the generated adsorption structures to the main database,
        and also to the corresponding database for the calculation type if specified

        Args:
        - adsorption_structures: list of AdsorptionStructure objects
        - filenames: list of strings with the filenames of the structures
        - calc_type (optional): if provided, the initial structures will be written
            also to the corresponding database for the calculation type
        '''

        # Write the adsorption structures to the database
        # excluding those that are already present
        calc_ids = []
        with ase.db.connect('structures.db') as db:
            for ads_struct, filename in zip(adsorption_structures, filenames):
                already_present = False
                for row in db.select():
                    if ads_struct.atoms == row.toatoms():
                        already_present = True
                        calc_ids.append(row.id)
                        break
                if not already_present:
                    calc_id = db.write(ads_struct.atoms,
                                       filename=filename,
                                       data=ads_struct.additional_data_arrays(),
                                       **ads_struct.to_info_dict())
                    calc_ids.append(calc_id)

        # Write the adsorption structures to the corresponding database
        if calc_type:
            with ase.db.connect(Database.calc_types[calc_type]) as db:
                for calc_id, ads_struct, filename in zip(calc_ids,adsorption_structures,filenames):
                    already_present = False
                    for row in db.select():
                        if ads_struct.atoms == row.toatoms():
                            break
                    if not already_present:
                        db.write(ads_struct.atoms,
                                calc_id=calc_id,
                                filename=filename,
                                **ads_struct.to_info_dict())


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


