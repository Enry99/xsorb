from ase.io import read, write
import numpy as np
import sys
#
from slab import Slab, adsorb_both_surfaces
from molecule import Molecule
from espresso_mod import Espresso_mod
from io_utils import get_energies, get_z, launch_jobs, restart_jobs
from settings import Settings
from filenames import *


def generate(SCF_RUN : bool, SAVEFIG=False, saveas_format=None): 
 
    #BEGIN STRUCTURES GENERATION ############################################################################

    settings=Settings()

    #Slab import from file
    slab = Slab(settings.slab_filename, layers_threshold=settings.layers_height, surface_sites_height=settings.surface_height)

    #sga = SpacegroupAnalyzer(slab.slab_pymat)
    #sops = sga.get_point_group_operations(cartesian=True)
    #for sop in sops:
    #    if sop.rotation_matrix[2][2] == 1.:
    #        print(sop.rotation_matrix)

    #Molecule import from file
    mol = Molecule(settings.molecule_filename, settings.molecule_axis_atoms, settings.axis_vector, settings.mol_subset_atoms)

    #Find adsorption sites and labels (site type and x,y coords.)
    adsites, adsites_labels = slab.find_adsorption_sites(
        *settings.sites_find_args.values(), 
        save_image=SAVEFIG,
        selected_sites=settings.selected_sites)
    for i, ads in enumerate(adsites):
        adsites[i][2] = max(slab.slab_ase.positions[:,2]) #takes care if some surface atoms are higher than the adsorption site

    distances_from_surf = len(settings.x_rot_angles)*[settings.scf_atom_distance]

    #Generate all the configs for the various molecular rotations and a list of labels
    all_mol_configs_ase, configs_labels = mol.generate_molecule_rotations(
        atom_index=settings.selected_atom_index,  
        vert_rotations=settings.y_rot_angles,
        screw_rotations=settings.x_rot_angles, 
        horiz_rotations=settings.z_rot_angles, 
        no_vert_rotx=settings.no_rot_vert,
        distance_from_surf=distances_from_surf, 
        min_distance=settings.scf_min_distance, 
        save_image=SAVEFIG
        )


    #Adsorption of molecule on all adsorption sites for all molecule orientations
    print('Generating adsorption structures...') 
    all_mol_on_slab_configs_ase = [slab.generate_adsorption_structures(molecule=mol_config, adsites=adsites) for mol_config in all_mol_configs_ase]
    all_mol_on_slab_configs_ase = sum(all_mol_on_slab_configs_ase, []) #flatten the 2D array
    if(False): all_mol_on_slab_configs_ase = adsorb_both_surfaces(all_mol_on_slab_configs_ase) #Replicate molecule on the other side of the slab. NOTE: Currently not working!
    full_labels = [mol_config[0]+site_label+mol_config[1] for mol_config in configs_labels for site_label in adsites_labels]
    print('All slab+adsorbate cells generated.')
       
    print('Writing pwi(s)...')
    csvfile=open(scf_labels_filename, 'w')
    csvfile.write('Label' + ',' + 'xrot' + ',' + 'yrot' + ',' + 'zrot' + ',' + 'site' + ',' + 'x' + ',' + 'y' + ',' + 'z' +'\n')

    if SCF_RUN or 'calculation' not in settings.espresso_settings_dict['CONTROL']: 
        settings.espresso_settings_dict['CONTROL'].update({'calculation' : 'scf'})
    settings.espresso_settings_dict['CONTROL'].update({'restart_mode' : 'from_scratch'})

    
    pwi_names = []

    for i in np.arange(len(all_mol_on_slab_configs_ase)):
        filename = pwi_prefix+'scf_'+str(i)+'.pwi'
        pwi_names.append(filename)
        calc = Espresso_mod(pseudopotentials=settings.pseudopotentials, 
                    input_data=settings.espresso_settings_dict,
                    filename=filename,
                    kpts= settings.kpoints[1] if 'gamma' not in settings.kpoints[0] else None, koffset=settings.kpoints[2] if 'gamma' not in settings.kpoints[0] else None)
        if(settings.fixed_layers_slab): 
            fixed_slab = slab.get_atoms_by_layers(settings.fixed_layers_slab)
        else: fixed_slab = settings.fixed_indices_slab
        calc.set_fixed_atoms(fixed_slab, slab.reindex_map, settings.fixed_indices_mol, mol.reindex_map, slab.natoms, mol.natoms, settings.fix_slab_xyz, settings.fix_mol_xyz)
        calc.set_system_flags(settings.starting_mag, settings.flags_i)
        calc.write_input(all_mol_on_slab_configs_ase[i])
        if(saveas_format is not None): write(filename.split('.')[0]+'.'+saveas_format, all_mol_on_slab_configs_ase[i])

        csvfile.write(str(i)+','+full_labels[i]+'\n')

    csvfile.close()

    print('All pwi(s) written.')
    #END OF STRUCTURE GENERATIONS #########################################################################################
    
    if SCF_RUN:
        launch_jobs(jobscript=settings.jobscript, pwi_list=pwi_names, outdirs=scf_outdir, jobname_prefix='scf', pwi_prefix=pwi_prefix, pwo_prefix=pwo_prefix)


def final_relax(threshold : float = None, exclude : list[int] = None, indices : list[int] = None):

    settings=Settings()
    

    #Slab import from file
    slab = Slab(settings.slab_filename, layers_threshold=settings.layers_height, surface_sites_height=settings.surface_height)

    #Molecule import from file
    mol = Molecule(settings.molecule_filename, settings.molecule_axis_atoms, settings.axis_vector, settings.mol_subset_atoms)

    if indices is not None: calcs = indices
    else:  
        print('Collecting energies from scf screening...')
        energies = get_energies(scf_labels_filename, scf_energies_filename, E_slab_mol=settings.E_slab_mol, pwo_prefix=pwo_prefix+'scf')
        if None in energies:
            print('Not all the calculations have reached convergence: impossible to identify the minimum. Quitting.')
            sys.exit(1)
        print('Scf energies collected.')


        if exclude is None: exclude = []
        else: print('Configurations {0} will be excluded, as requested'.format(exclude))

        e_min = min([energies[i] for i in [*range(len(energies))] if i not in exclude])
        i_minimum = energies.index(e_min)
        calcs = []
        subset_energies = []
        if threshold is not None:
            for i, energy in enumerate(energies):
                if i in exclude: continue
                if energy - e_min <= threshold:
                    calcs.append(i)
                    subset_energies.append(energy)
            #sort them by minimum energy, in order to launch first the minimum
            sorted_indices = np.argsort(subset_energies, kind='stable').tolist()
            calcs = [calcs[i] for i in sorted_indices]
        else: calcs = [i_minimum]


    settings.espresso_settings_dict['CONTROL'].update({'calculation' : 'relax'})
    settings.espresso_settings_dict['CONTROL'].update({'restart_mode' : 'from_scratch'})
    settings.espresso_settings_dict['IONS'].update({'ion_dynamics': settings.ion_dynamics})


    #Re-generation of structures, placing the molecule very close to the surface to avoid
    #trapping in local physisorption minima
    adsites, adsites_labels = slab.find_adsorption_sites(
        *settings.sites_find_args.values(), 
        save_image=False,
        selected_sites=settings.selected_sites)

    all_mol_configs_ase, configs_labels = mol.generate_molecule_rotations(
        atom_index=settings.selected_atom_index,  
        vert_rotations=settings.y_rot_angles,
        screw_rotations=settings.x_rot_angles, 
        horiz_rotations=settings.z_rot_angles, 
        no_vert_rotx=settings.no_rot_vert,
        distance_from_surf=settings.rel_atom_distance, 
        min_distance=settings.rel_min_distance, 
        save_image=False
        )

    #Adsorption of molecule on all adsorption sites for all molecule orientations
    print('Generating adsorption structures...') 
    all_mol_on_slab_configs_ase = [slab.generate_adsorption_structures(molecule=mol_config, adsites=adsites) for mol_config in all_mol_configs_ase]
    all_mol_on_slab_configs_ase = sum(all_mol_on_slab_configs_ase, []) #flatten the 2D array
    if(False): all_mol_on_slab_configs_ase = adsorb_both_surfaces(all_mol_on_slab_configs_ase) #Replicate molecule on the other side of the slab. NOTE: Currently not working!
    full_labels = [mol_config[0]+site_label+mol_config[1] for mol_config in configs_labels for site_label in adsites_labels]
    print('All slab+adsorbate cells generated.')

    pwi_names = []
    for i in calcs:       
        #struct_ase = read(pwi_prefix+'scf_'+str(i)+'.pwi') #simply reads the files, avoid to re-generate them

        fixed_indices_molecule = slab.natoms + np.array(settings.fixed_indices_mol)

        filename = pwi_prefix+'relax_'+str(i)+'.pwi'
        pwi_names.append(filename)
        calc = Espresso_mod(pseudopotentials=settings.pseudopotentials, 
                    input_data=settings.espresso_settings_dict,
                    filename=filename,
                    kpts= settings.kpoints[1] if 'gamma' not in settings.kpoints[0] else None, koffset=settings.kpoints[2] if 'gamma' not in settings.kpoints[0] else None)
        if(settings.fixed_layers_slab): fixed_slab = slab.get_atoms_by_layers(settings.fixed_layers_slab)
        else: fixed_slab = settings.fixed_indices_slab
        calc.set_fixed_atoms(
            fixed_slab,
            slab.reindex_map, 
            fixed_indices_molecule,
            mol.reindex_map, 
            slab.natoms, 
            mol.natoms,
            settings.fix_slab_xyz, settings.fix_mol_xyz)
        calc.set_system_flags(settings.starting_mag, settings.flags_i)
        calc.write_input(all_mol_on_slab_configs_ase[i])

    launch_jobs(jobscript=settings.jobscript, pwi_list=pwi_names, outdirs=relax_outdir, jobname_prefix='rel', pwi_prefix=pwi_prefix, pwo_prefix=pwo_prefix)