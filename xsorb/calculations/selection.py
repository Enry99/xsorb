#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Enrico Pedretti

'''
Module with the functions to select the configurations for subsequent calculations
(both from the screening and from the ML optimization)
'''

from __future__ import annotations
import logging

from xsorb.ase_custom.atoms import AtomsCustom
from xsorb.io.database import Database
from xsorb.adsorptiondata.adsorptionstructure import AdsorptionStructure


def _select_calculations(rows : list,
                        n_configs : int | None = None,
                        threshold : float | None = None) -> list:
    '''
    Returns the indices of the configurations to be relaxed, according to the specified criteria.

    Args:
    - rows: list of rows from the database, already sorted by energy
    - n_configs: number of configurations to be relaxed, starting from the one with lowest energy
    - threshold: energy threshold (in eV) from the NOT EXCLUDED lowest energy configuration.
        The configuration with E - Emin < threshold will be selected

    Returns:
    - calculations_indices: list of indices of the configurations to be relaxed
    '''

    if threshold is not None:
        emin = [row.energy for row in rows][0]
        calc_indices = [row.get('calc_id') for row in rows if row.energy - emin < threshold]

    elif n_configs is not None:
        calc_indices = [row.get('calc_id') for i, row in enumerate(rows) if i < n_configs]

    else:
        raise ValueError('Either n_configs or threshold must be specified.')

    return calc_indices


def obtain_calc_indices(*,
                        calc_type : str,
                        n_configs: int | None = None,
                        threshold : float | None = None,
                        excluded_calc_ids : list | None = None,
                        by_site : bool = False,
                        by_mol_atom : bool = False,
                        separate_chem_phys : bool = False,
                        verbose : bool = True) -> list[int]:
    '''
    Returns a list with the indices of the configurations from previous calculations,
    according to the specified criteria.
    If no criteria is specified, all the configurations are returned.

    Args:
    - calc_type: type of calculation to get the indices from. Can be 'screening', 'mlopt'
    - n_configs: number of configurations to be relaxed, starting from the one with lowest energy
    - threshold: energy threshold (in eV) from the NOT EXCLUDED lowest energy configuration.
        The configuration with E - Emin < threshold will be selected
    - excluded_calc_ids: indices of the configurations to be excluded
    - by_site: do the configuration identification separately for each site.
    - by_mol_atom: do the configuration identification separately for each ref. atom of the molecule
    - separate_chem_phys: do the configuration identification separately
        for physisorption and chemisorption
    - verbose: print messages

    Returns:
    - selected_calc_ids: list of indices of the configurations to be relaxed
    '''

    # check the input parameters
    if n_configs is not None and threshold is not None:
        raise ValueError('Only one between n_configs and threshold can be specified.')

    if verbose:
        logging.info(f'Collecting results from {calc_type}...')
        if excluded_calc_ids is not None:
            logging.info(f'Configurations {excluded_calc_ids} will be excluded, as requested.')


    #### start collecting the results from the database ####

    # case 1: no criteria specified, just return all the indices of the configurations
    if n_configs is None and threshold is None:
        rows = Database.get_calculations(calc_type=calc_type,
                                        exclude_ids=excluded_calc_ids,
                                        include_data=False)
        return [row.get('calc_id') for row in rows]


    # case 2: some criteria specified, so we need to select the configurations

    #select only rows that have energy, sorting them from lowest to highest energy
    rows = Database.get_calculations(calc_type=calc_type,
                                    selection='energy',
                                    exclude_ids=excluded_calc_ids,
                                    #columns=['energy', 'calc_id', 'site', 'bonds'],
                                    sort_key='energy',
                                    include_data=False)


    # helper functions to make the code a bit cleaner
    def _chemphys_subsets(rows_list: list, choose_by_subsets: bool) -> list[list]:
        '''
        Returns either [chemisorption_rows, physisorption_rows]
        or [all_rows] depending on the choose_by_subsets flag.
        '''
        if choose_by_subsets:
            return [[row for row in rows_list if row.bonds != 'none'],
                        [row for row in rows_list if row.bonds == 'none']]
        else:
            return [rows_list]

    def _mol_atoms_subsets(rows_list: list, choose_by_subsets: bool) -> list:
        '''
        Returns either [rows_mol_atom1, rows_mol_atom2, ...]
        or [all_rows] depending on the choose_by_subsets flag.
        '''
        if choose_by_subsets:
            local_rows = []
            for mol_atom in set(row.get('mol_atom') for row in rows_list):
                local_rows.append([row for row in rows_list if row.get('mol_atom') == mol_atom])
            return local_rows
        else:
            return [rows_list]


    def _site_subsets(rows_list: list, choose_by_subsets: bool) -> list:
        '''
        Returns either [rows_site1, rows_site2, ...]
        or [all_rows] depending on the choose_by_subsets flag.
        '''
        if choose_by_subsets:
            local_rows = []
            for site in set(row.get('site') for row in rows_list):
                local_rows.append([row for row in rows_list if row.get('site') == site])
            return local_rows
        else:
            return [rows_list]


    # loop over the subsets, adding the ids as a flat list
    selected_calc_ids = []
    for subset1 in _mol_atoms_subsets(rows, by_mol_atom):
        for subset2 in _site_subsets(subset1, by_site):
            for subset3 in _chemphys_subsets(subset2, separate_chem_phys):

                calc_ids = _select_calculations(rows=subset3,
                                                n_configs=n_configs,
                                                threshold=threshold)

                selected_calc_ids.extend(calc_ids)

    # DEBUG: check that no duplicates are present
    if len(selected_calc_ids) != len(set(selected_calc_ids)):
        raise ValueError('Duplicate calculation IDs found in the selected configurations.')


    if verbose: logging.info(f'{calc_type} results collected.') #pylint: disable=multiple-statements

    return selected_calc_ids


def get_adsorption_structures(get_structures_from : str,
                              calc_ids : list[int] | None = None) -> list[AdsorptionStructure]:
    '''
    Returns the a list of AdsorptionStructure objects, but with the atoms
    substituted with the ones from the previous calculation.
    If from mlopt, the the original constraints for dft are set, important when
    the mlopt was performed by fixing the slab.

    Args:
    - calc_ids: list of indices of the configurations to be retrieved. If None, all are retrieved
    - get_structures_from: type of calculation to get the structures from.
        Can be 'screening', 'mlopt', 'structures'

    Returns:
    - adsorption_structures: list of AdsorptionStructure objects
    '''

    #get structures from database
    rows = Database.get_calculations(calc_type=get_structures_from, calc_ids=calc_ids)
    rows_original = Database.get_structures(calc_ids=calc_ids)


    # prepare the structures by substituting the atoms with the one from the previous calculation
    adsorption_structures : list [AdsorptionStructure] = []

    for row, row_original in zip(rows, rows_original):
        ads_struct = AdsorptionStructure.fromdict(
            row.data.AdsorptionCalculation.get('adsorption_structure'))
        atoms = AtomsCustom(row.toatoms())

        if get_structures_from == 'mlopt':
            atoms.set_constraint(row_original.get('constraints'))

        ads_struct.atoms = atoms
        adsorption_structures.append(ads_struct)

    return adsorption_structures
