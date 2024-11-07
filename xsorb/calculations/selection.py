#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Enrico Pedretti

'''
Module with the functions to select the configurations for subsequent calculations
(both from the screening and from the ML optimization)
'''

from __future__ import annotations
from typing import TYPE_CHECKING

from ase import Atoms

from xsorb.io.database import Database
if TYPE_CHECKING:
    from xsorb.structures.properties import AdsorptionStructure


def select_calculations(rows : list,
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
        calc_indices = [row.calc_id for row in rows if row.energy - emin < threshold]

    elif n_configs is not None:
        calc_indices = [row.calc_id for i, row in enumerate(rows) if i < n_configs]

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
    Returns a list with the indices of the configurations to be relaxed,
    according to the specified criteria.
    If no criteria is specified, all the configurations are returned.

    Args:
    - calc_type: type of calculation to get the indices from. Can be 'screening', 'ml_opt'
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

    if n_configs is None and threshold is None:
        #just return the indices of all the configurations
        rows, _ = Database.get_calculations(calc_type=calc_type,
                                        exclude_ids=excluded_calc_ids,
                                        columns=['calc_id'],
                                        include_data=False)
        return [row.calc_id for row in rows]

    if n_configs is not None and threshold is not None:
        raise ValueError('Only one between n_configs and threshold can be specified.')

    if verbose: print(f'Collecting results from {calc_type}...') #pylint: disable=multiple-statements

    if excluded_calc_ids is not None and verbose:
        print(f'Configurations {excluded_calc_ids} will be excluded, as requested.')


    selected_calc_ids = []
    #select only rows that have energy
    selections = ['energy,bonds!=None', 'energy,bonds=None'] if separate_chem_phys else ['energy']
    for selection in selections:
        rows, _ = Database.get_calculations(calc_type=calc_type,
                                        selection=selection,
                                        exclude_ids=excluded_calc_ids,
                                        columns=['energy', 'calc_id', 'site'],
                                        sort_key='energy',
                                        include_data=False)

        if by_mol_atom:
            for mol_atom in set(row.mol_atom for row in rows):
                rows_mol_atom = [row for row in rows if row.mol_atom  == mol_atom]

                if by_site:
                    for site in set(row.site for row in rows_mol_atom):
                        rows_site = [row for row in rows_mol_atom if row.site == site]
                        selected_calc_ids += select_calculations(rows_site, n_configs, threshold)
                else:
                    selected_calc_ids += select_calculations(rows_mol_atom, n_configs, threshold)
        else:
            if by_site:
                for site in set(row.site for row in rows):
                    rows_site = [row for row in rows if row.site == site]
                    selected_calc_ids += select_calculations(rows_site, n_configs, threshold)
            else:
                selected_calc_ids += select_calculations(rows, n_configs, threshold)

    print(f'{calc_type} results collected.')

    return selected_calc_ids


def get_adsorption_structures(get_structures_from : str,
                              calc_ids : list[int] | None = None) -> list[AdsorptionStructure]:
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
    rows, _ = Database.get_calculations(calc_type=get_structures_from, calc_ids=calc_ids)

    if get_structures_from == 'ml_opt':
        rows_original, _ = Database.get_calculations(calc_type='structures',
                                                  calc_ids=calc_ids)
        constraints = [row.constraints for row in rows_original]

    # prepare the structures by substituting the atoms with the one from the
    # previous calculation
    adsorption_structures : list [AdsorptionStructure] = []
    for row in rows:
        ads_struct = AdsorptionStructure(**row.data.adsorption_structure)
        atoms : Atoms = row.toatoms()

        if get_structures_from == 'ml_opt':
            atoms.set_constraint(constraints.pop(0))

        ads_struct.atoms = atoms
        adsorption_structures.append(ads_struct)

    return adsorption_structures
