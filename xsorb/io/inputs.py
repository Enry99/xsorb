'''
Module for writing input files for the calculations
'''

from xsorb.settings import Settings
from xsorb.structures import AdsorptionStructure
from xsorb.io.database import Database


#WORKFLOW:
#per screening from_scratch, dobbiamo generare nuove strutture.
#queste strutture non hanno ancora indici, quindi dobbiamo scorrere il
#database delle strutture per vedere quali indici sono già stati assegnati,
#usare quelli e assegnare gli altri.

#gestire caso dei file già








def write_inputs(settings : Settings,
                 calc_type : str,
                 adsorption_structures : list[AdsorptionStructure],
                 calc_indices : list[int] | None = None,
                 override_settings : bool = True,
                 ask_before_overwrite : bool = False,
                 verbose : bool = True) -> list[int]:
    '''
    Writes the input files for all the adsorption configurations.


    Args:
    - settings: Settings object, containing all the parameters
    - adsorption_structures: list AdsorptionStructure objects
    - calc_type: 'screening', 'relax' or 'ml_opt'
    - override_settings: override some specifc settings (e.g. conv tresholds)
    - interactive: interactive mode: ask before overwriting files that are already present

    Returns:
    - written_indices: list of the indices of the configurations
        for which the input files have been written
    '''

    if verbose: print('Writing input files...') #pylint: disable=multiple-statements

    if calc_indices is None:

    #when launching ml or screening from scratch, write to db_initial.
    #when launching screenining from preopt, read from db_ml, do not write to db_initial
    #when launching relax, read either from db_ml or db_screening.



    # ANSWER_ALL = False
    # answer = 'yes'

    written_indices = []


    if override_settings: override_settings(settings, calc_type)

    for i, atoms in zip(calc_indices, all_mol_on_slab_configs_ase):

        file_label = f'{calc_type.lower()}_{i}'

        corresponding_outfile = OUT_FILE_PATHS[calc_type][settings.program].format(i)

        #TODO: for now, always overwrite. Implement interactive mode later
        #if in interactive mode, ask before overwriting
        # if(interactive and os.path.isfile(corresponding_outfile)): #convoluted, but it works
        #     print(f'{corresponding_outfile} already present, possibly from a running calculation. '+ \
        #           ('It will be {0}, as requested.'.format('skipped' if 'n' in answer  else 're-calculated') \
        #           if ANSWER_ALL else 'You can decide to re-calculate it or skip it.'))
        #     while True and not ANSWER_ALL:
        #         answer = input('Re-calculate? ("y" = yes to this one, "yall" = yes to all, "n" = no to this one, "nall" = no to all): ')
        #         if answer == 'yes' or answer == 'y' or answer == 'yall' or answer == 'no' or answer == 'n' or answer == 'nall':
        #             if answer == 'yall' or answer == 'nall': ANSWER_ALL = True
        #             break
        #         else: print('Value not recognized. Try again.')
        #     if answer == 'no' or answer == 'n' or answer == 'nall': continue #skip if user does not want to overwrite


        j_dir = f'{screening_outdir if calc_type == "SCREENING" else relax_outdir}/{i}'
        calc = Calculator(settings, file_label, atoms, j_dir)
        calc.write_input(atoms)
        written_indices.append(i)

    if interactive: print('All input files written.')

    return written_indices


def write_inputs_ml(settings : Settings,
                 all_mol_on_slab_configs_ase : list,
                 calc_indices : list,
                 VERBOSE : bool = True):
    '''
    Writes the input files for all the adsorption configurations.


    Args:
    - settings: Settings object, containing the dft code parameters
    - all_mol_on_slab_configs_ase: list of ASE atoms for all the adsorption configurations

    Returns:
    - written_indices: list of the indices of the configurations for which the input files have been written
    '''

    if VERBOSE: print('Writing pre-optimization input files...')

    written_indices = []

    os.makedirs(preopt_outdir, exist_ok=True)

    for i, atoms in zip(calc_indices, all_mol_on_slab_configs_ase):

        os.makedirs(f'{preopt_outdir}/{i}', exist_ok=True)
        atoms.pbc = True #ensure that the periodic boundary conditions are set
        write(IN_FILE_PATHS['PREOPT']['ML'].format(i), atoms)
        written_indices.append(i)

    if VERBOSE: print('All input files written.')

    return written_indices