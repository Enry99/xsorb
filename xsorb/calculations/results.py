'''
Module to check and update the results of the calculations
'''

from dataclasses import dataclass
from typing import Optional
from pathlib import Path

from ase import Atoms

from xsorb.structures.utils import slab_mol_bonds
from xsorb.io.utils import ase_custom_read as read
from xsorb.io.settings import Settings
from xsorb.dft_codes.definitions import OUT_FILE_PATHS, LOG_FILE_PATHS, \
    SCF_NONCONVERGED_STRINGS, SCF_CONVERGED_STRINGS, OPTIMIZATION_COMPLETED_STRINGS
from xsorb.io.database import Calc_Paths


@dataclass
class CalculationResults:
    '''
    Dataclass to store the results of a calculation
    '''
    atoms: Atoms
    adsorption_energy: float
    status : str #'completed', 'incomplete'
    scf_nonconverged : bool
    bonds : str
    trajectory : Optional[list[Atoms]]

    @property
    def energy(self):
        '''
        Obtain the energy of the atoms object
        '''
        return self.atoms.get_potential_energy()

    @property
    def energy_evolution(self):
        '''
        Obtain the energy evolution along the trajectory
        '''
        if self.trajectory is not None:
            return [atoms.get_potential_energy() for atoms in self.trajectory]
        else:
            return None


def read_total_e_slab_mol_from_default_folders(ml: bool, verbose : bool = True):
    '''
    Attempt to read the energies of the slab and molecule from the
    output files of the slab/molecule calculations launched by xsorb
    (not provided by the user in the settings file)

    Args:
    - ml: if True, the calculation is for the ML optimization
    - verbose: print messages

    '''
    try:
        slab_en = read(self.input.slab_filename).get_potential_energy()
        mol_en = read(self.input.molecule_filename).get_potential_energy()
        self.input.E_slab_mol = [slab_en, mol_en]
    except Exception as e: # pylint: disable=broad-except
        #this is a general exception to catch any error that might occur.
        #It can be file not found, or energy not present in the file, etc.
        self.input.E_slab_mol = [0.0, 0.0]
        if verbose:
            print('It was not possible to obtain slab and molecule energy in any way.',
                    f'Error message from ase: {e}.',
                    'Total energies will be shown instead of adsorption energies.')


def is_optimization_completed(filename : str, program : str):
    '''
    Check if the given calculation is completed, reading the output file

    Args:
    - filename: path to the LOG_FILE(==output file for espresso)
    - calc_type: 'screening','relax','ml_opt'

    Returns:
    True or False
    '''

    searchfor = OPTIMIZATION_COMPLETED_STRINGS[program]

    with open(filename, 'r') as f:
        file_content = f.readlines()

    completed = False
    for line in file_content:
        if searchfor in line:
            completed = True
            break

    return completed


def is_scf_not_converged(filename : str, program : str):
    '''
    Check if the given calculation is completed, reading the output file

    Args:
    - filename: path to the LOG_FILE(==output file for espresso)
    - calc_type: 'screening','relax','ml_opt'

    Returns:
    True or False
    '''

    if program == 'ml': return False #pylint: disable=multiple-statements

    searchfor = SCF_NONCONVERGED_STRINGS[program]
    convergence_string = SCF_CONVERGED_STRINGS[program]

    with open(filename, 'r') as f:
        file_content = f.readlines()

    # we might encounter the situation where a first loop is not converged,
    # but the last one is, so we need to check all the lines:
    # the last one (conv or not conv) determines the status
    nonconv = False
    for line in file_content:
        if searchfor in line:
            nonconv = True
        elif convergence_string in line:
            nonconv = False

    return nonconv


def get_atoms_from_calc(filename : str, return_trajectory : bool = True):
    '''
    Reads the output file and returns the atoms object.
    If there is an error reading the file, it prints a message and returns None

    Args:
    - filename: path to the output file
    - return_trajectory: if True, returns a list of Atoms objects

    Returns:
    The atom object from the output file, or a list of Atoms objects if return_trajectory is True
    '''

    if return_trajectory:
        try:
            trajectory = read(filename, index=':')
        except Exception as exc:
            print(f'Error reading trajectory from file {filename}: {exc}. '\
                'Attempting to read only the last configuration.')
            try:
                trajectory = [read(filename)]
            except Exception as exc2:
                print(f'Fatal error reading file {filename}: {exc2}.')
                return None

        return trajectory

    else:
        try:
            atoms = read(filename)
        except Exception as exc:
            print(f'Error reading file {filename}: {exc}.')
            return None

        return atoms


def get_calculations_results(paths_list: list[Calc_Paths],
                             program : str,
                             calc_type : str,
                             total_e_slab_mol : float | None,
                             mult : float,
                             verbose : bool =True):
    '''
    Reads the output files and returns a list of CalculationResults objects.
    The calculations with no output will be added as None

    Args:
    - paths_list: list of Paths objects, containg in_file_path, out_file_path and log_file_path
    - program: DFT program. Possible values: 'espresso','vasp','ml'
    - calc_type: 'screening', 'relax' or 'ml_opt'
    - total_e_slab_mol: total energy of the slab and molecule, if available
    - mult: multiplicative factor for the covalent radii to determine bonding.
    '''

    if total_e_slab_mol is None:
        total_e_slab_mol = read_total_e_slab_mol_from_calc()

    results_list : list[CalculationResults | None ] = []

    for paths in paths_list:
        if not Path(paths.out_file_path).exists() or not Path(paths.log_file_path).exists():
            if verbose:
                print(f'Warning! File {paths.out_file_path} not found. Skipping.')
            results_list.append(None)
            continue

        if is_optimization_completed(paths.log_file_path, program):
            status = 'completed'
        else:
            status = 'incomplete'

        scf_nonconverged = is_scf_not_converged(paths.log_file_path, program)
        if verbose and scf_nonconverged:
            print(f'Warning! {paths.out_file_path} failed to reach SCF convergence. '\
                  'Number of electronic steps exceeded.')

        traj = get_atoms_from_calc(paths.out_file_path)

        if traj:
            atoms = traj[-1]
            calc_results = CalculationResults(atoms=atoms,
                                              status=status,
                                              scf_nonconverged=scf_nonconverged,
                                              bonds=None,
                                              trajectory=traj)
        else:
            calc_results = None

        results_list.append(calc_results)

    return results_list




def write_results_to_file(TXT=False):
    '''
    Function to write the calculations results to a csv file.
    It can be called before all the jobs have finished.

    It needs to read the 'site_labels.csv' file.

    Args:
    - TXT: write a txt file (tab separated) instead of csv, sorted by screening or relax energies
    '''
    #

    settings = Settings()

    datafile = pd.read_csv(labels_filename, index_col=0)


    if os.path.isdir(preopt_outdir): #for preopt
            preopt_results = \
                get_calculations_results_ml()

            column_name = 'Eads_pre(eV)'

            column_data = []
            for i in datafile.index: #if the file exists, and so the energy might be present, or it might be None
                if i in preopt_results['energies']:
                    if preopt_results['energies'][i] is None:
                        print(f'Config. {i} has no energy. It will be skipped.')
                        column_data.append(None)
                        continue

                    column_data.append(preopt_results['energies'][i] )

                    if not preopt_results['relax_completed'][i]:
                        print(f'Warning! {i} relaxation has not reached final configuration. The energy will be marked with a *')
                        column_data[-1] = f'{column_data[-1]:.3f}*'

                else: #if the file does not exist
                    column_data.append(None)

            datafile[column_name] = column_data


    if os.path.isdir(screening_outdir): #for screening
        screening_results = \
            get_calculations_results(program=settings.program, calc_type='SCREENING', E_slab_mol=settings.E_slab_mol, VERBOSE=False)

        column_name = 'Eads_scr(eV)' if 0 not in settings.E_slab_mol else 'Etot_scr(eV)'

        column_data = []
        for i in datafile.index: #if the file exists, and so the energy might be present, or it might be None
            if i in screening_results['energies']:
                if screening_results['energies'][i] is None:
                    print(f'Config. {i} has not reached the first scf convergence. It will be skipped.')
                    column_data.append(None)
                    continue

                column_data.append(screening_results['energies'][i] )

                if screening_results['scf_nonconverged'][i]:
                    print(f'Warning! {i} failed to reach SCF convergence in the last ionic step. The energy will be marked with **')
                    column_data[-1] = f'{column_data[-1]:.3f}**'
                elif not screening_results['relax_completed'][i]:
                    print(f'Warning! {i} relaxation has not reached final configuration. The energy will be marked with a *')
                    column_data[-1] = f'{column_data[-1]:.3f}*'

            else: #if the file does not exist
                column_data.append(None)

        datafile[column_name] = column_data


    if os.path.isdir(relax_outdir): #for relax
        relax_results = \
            get_calculations_results(program=settings.program, calc_type='RELAX', E_slab_mol=settings.E_slab_mol, VERBOSE=False)

        column_name = 'Eads_rel(eV)' if 0 not in settings.E_slab_mol else 'Etot_rel(eV)'

        column_data = []
        for i in datafile.index: #if the file exists, and so the energy might be present, or it might be None
            if i in relax_results['energies']:
                if relax_results['energies'][i] is None:
                    print(f'Config. {i} has not reached the first scf convergence. It will be skipped.')
                    column_data.append(None)
                    continue

                column_data.append(relax_results['energies'][i] )

                if relax_results['scf_nonconverged'][i]:
                    print(f'Warning! {i} failed to reach SCF convergence in the last ionic step. The energy will be marked with **')
                    column_data[-1] = f'{column_data[-1]:.3f}**'
                elif not relax_results['relax_completed'][i]:
                    print(f'Warning! {i} relaxation has not reached final configuration. The energy will be marked with a *')
                    column_data[-1] = f'{column_data[-1]:.3f}*'

            else: #if the file does not exist
                column_data.append(None)

        datafile[column_name] = column_data

    #add bonding info
    slab = Slab(settings.slab_filename)
    mol = Molecule(settings.molecule_filename, atom_index=settings.selected_atom_index)
    mol_indices = np.arange(mol.natoms) if settings.mol_before_slab else np.arange(mol.natoms) + slab.natoms
    bonding_status = []
    for i in datafile.index: #if the file exists, and so the energy might be present, or it might be None
        if os.path.isdir(relax_outdir) and i in relax_results['energies'] and relax_results['energies'][i]: #prioritize status from relax over screening
            status = check_bond_status(settings.program,
                                                    calc_type='RELAX',
                                                    i_calc=i,
                                                    mol_indices=mol_indices)
            bonding_status.append(','.join(status) if status else 'No')
        elif os.path.isdir(screening_outdir) and i in screening_results['energies'] and screening_results['energies'][i]:
            status = check_bond_status(settings.program,
                                                    calc_type='SCREENING',
                                                    i_calc=i,
                                                    mol_indices=mol_indices)
            bonding_status.append(','.join(status) if status else 'No')
        elif os.path.isdir(preopt_outdir) and i in preopt_results['energies'] and preopt_results['energies'][i]:
            status = check_bond_status('ML',
                                                    calc_type='PREOPT',
                                                    i_calc=i,
                                                    mol_indices=mol_indices)
            bonding_status.append(','.join(status) if status else 'No')
        else: #if the file does not exist
            bonding_status.append(None)

    datafile['Bonded'] = bonding_status


    if(TXT): #sort by energy column (relax if available, else screening)
        datafile.sort_values(by=column_name)

    datafile.to_csv(results_filename.replace('csv', 'txt' if TXT else 'csv'), sep='\t' if TXT else ',')

    print('Results file written.')


def check_bond_status(program : str, calc_type : str, i_calc : int, mol_indices : list):
    '''
    Reads output file and returns the list of bonds between slab and molecule

    Args:
    - program: DFT program. Possible values: 'ESPRESSO' or 'VASP'
    - calc_type: 'SCREENING' or 'RELAX'
    - i_calc: numeric index of the calculation
    - mol_indices: indices of the atoms belonging to the molecule
    '''
    filename = OUT_FILE_PATHS[calc_type][program].format(i_calc)
    atoms = read(filename)

    slab = atoms[[atom.index for atom in atoms if atom.index not in mol_indices]]
    mol = atoms[mol_indices]

    return slab_mol_bonds(slab, mol)