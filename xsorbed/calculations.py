from ase.io import read, write
import numpy as np
import os, sys
from natsort import natsorted
import glob
#
from slab import Slab, adsorb_both_surfaces
from molecule import Molecule
from espresso_mod import Espresso_mod
from io_utils import get_energies, launch_jobs
from settings import Settings
from filenames import *


def generate(RUN : bool, etot_forc_conv = [5e-3, 5e-2], SAVEFIG=False, saveas_format=None): 
 
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


    #Generate all the configs for the various molecular rotations and a list of labels
    all_mol_configs_ase, configs_labels = mol.generate_molecule_rotations(
        atom_index=settings.selected_atom_index,  
        y_rot_angles=settings.y_rot_angles,
        x_rot_angles=settings.x_rot_angles, 
        z_rot_angles=settings.z_rot_angles, 
        vert_angles_list=settings.vertical_angles,
        distance_from_surf=settings.screening_atom_distance, 
        min_distance=settings.screening_min_distance, 
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
    csvfile=open(labels_filename, 'w')
    csvfile.write('Label' + ',' + 'xrot' + ',' + 'yrot' + ',' + 'zrot' + ',' + 'site' + ',' + 'x' + ',' + 'y' + ',' + 'z' +'\n')

    if RUN: 
        settings.espresso_settings_dict['CONTROL'].update({'calculation' : 'relax' })
        settings.espresso_settings_dict['CONTROL'].update({'etot_conv_thr' : etot_forc_conv[0]})
        settings.espresso_settings_dict['CONTROL'].update({'forc_conv_thr' : etot_forc_conv[1]})
        settings.espresso_settings_dict['IONS'].update({'upscale': 1})
    settings.espresso_settings_dict['CONTROL'].update({'restart_mode' : 'from_scratch'})

    if saveas_format is not None:
        folder = saveas_format+'/screening/'
        if not os.path.exists(saveas_format):
            os.mkdir(saveas_format)
        if not os.path.exists(folder):
            os.mkdir(folder)  

    
    pwi_names = []

    for i in np.arange(len(all_mol_on_slab_configs_ase)):
        filename = pw_files_prefix+'screening_'+str(i)+'.pwi'
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
        if(saveas_format is not None): write(folder+filename.split('.')[0]+'.'+saveas_format, all_mol_on_slab_configs_ase[i])

        csvfile.write(str(i)+','+full_labels[i]+'\n')

    csvfile.close()

    print('All pwi(s) written.')
    #END OF STRUCTURE GENERATIONS #########################################################################################
    
    if RUN:
        launch_jobs(jobscript=settings.jobscript, pwi_list=pwi_names, outdirs=screening_outdir, jobname_title='scr')


def final_relax(n_configs: int = None, threshold : float = None, exclude : list= None, indices : list = None, REGENERATE=False):
    
    if n_configs is None and threshold is None: n_configs = N_relax_default
    
    #Check for duplicates in exclude or in indices
    if indices is not None:
        if (len(set(indices)) < len(indices)):
            print('The indices list contain duplicate elements. Quitting.')
            sys.exit(1)
    if exclude is not None:
        if (len(set(exclude)) < len(exclude)):
            print('The exclude list contain duplicate elements. Quitting.')
            sys.exit(1)        


    settings=Settings()
    
    if indices is not None: calcs = indices
    else:
        print('Collecting energies from screening...')
        energies = get_energies(E_slab_mol=settings.E_slab_mol, pwo_prefix='screening')
        if None in energies: #TODO: sistemare!!!! per ora non differenzia lo screening ibrido e richiede solo che si sia raggiunta la conv scf
            print('Not all the calculations have reached convergence: impossible to identify the minimum. Quitting.')
            sys.exit(1)
        print('screening energies collected.')


        if exclude is None: exclude = []
        else: print('Configurations {0} will be excluded, as requested'.format(exclude))


        calcs = []
        subset_energies = []
        if n_configs is not None:
            sorted_indices = np.argsort(energies, kind='stable').tolist()
            #print(sorted_indices)
            for i in sorted_indices:
                if i in exclude: continue
                if len(calcs)<n_configs:
                    calcs.append(i)
        elif threshold is not None:
            e_min = min(energies)
            for i, energy in enumerate(energies):
                if i in exclude: continue
                if energy - e_min <= threshold:
                    calcs.append(i)
                    subset_energies.append(energy)
            #sort them by minimum energy, in order to launch first the minimum. TODO: do this also for the --n case
            sorted_indices = np.argsort(subset_energies, kind='stable').tolist()
            calcs = [calcs[i] for i in sorted_indices]


    settings.espresso_settings_dict['CONTROL'].update({'calculation' : 'relax'})
    settings.espresso_settings_dict['CONTROL'].update({'restart_mode' : 'from_scratch'})
    settings.espresso_settings_dict['IONS'].update({'ion_dynamics': settings.ion_dynamics})

    #Slab import from file
    slab = Slab(settings.slab_filename, layers_threshold=settings.layers_height, surface_sites_height=settings.surface_height)

    #Molecule import from file
    mol = Molecule(settings.molecule_filename, settings.molecule_axis_atoms, settings.axis_vector, settings.mol_subset_atoms)


    if(REGENERATE):
        #Re-generation of structures, placing the molecule very close to the surface to avoid
        #trapping in local physisorption minima
        adsites, adsites_labels = slab.find_adsorption_sites(
            *settings.sites_find_args.values(), 
            save_image=False,
            selected_sites=settings.selected_sites)

        all_mol_configs_ase, configs_labels = mol.generate_molecule_rotations(
            atom_index=settings.selected_atom_index,  
            y_rot_angles=settings.y_rot_angles,
            x_rot_angles=settings.x_rot_angles, 
            z_rot_angles=settings.z_rot_angles, 
            vert_angles_list=settings.vertical_angles,
            distance_from_surf=settings.relax_atom_distance, 
            min_distance=settings.relax_min_distance, 
            save_image=False
            )

        #Adsorption of molecule on all adsorption sites for all molecule orientations
        print('Generating adsorption structures...') 
        all_mol_on_slab_configs_ase = [slab.generate_adsorption_structures(molecule=mol_config, adsites=adsites) for mol_config in all_mol_configs_ase]
        all_mol_on_slab_configs_ase = sum(all_mol_on_slab_configs_ase, []) #flatten the 2D array
        if(False): all_mol_on_slab_configs_ase = adsorb_both_surfaces(all_mol_on_slab_configs_ase) #Replicate molecule on the other side of the slab. NOTE: Currently not working!
        full_labels = [mol_config[0]+site_label+mol_config[1] for mol_config in configs_labels for site_label in adsites_labels]
        print('All slab+adsorbate cells generated.')
    else:
        files = natsorted(glob.glob( pw_files_prefix + "screening_*.pwo" ))        
        all_mol_on_slab_configs_ase = [read(filename=file, results_required=False) for file in files]

 

    pwi_names = []
    for i in calcs:       
        #struct_ase = read(pwi_prefix+'screening_'+str(i)+'.pwi') #simply reads the files, avoid to re-generate them

        fixed_indices_molecule = slab.natoms + np.array(settings.fixed_indices_mol)

        filename = pw_files_prefix+'relax_'+str(i)+'.pwi'
        pwi_names.append(filename)
        calc = Espresso_mod(pseudopotentials=settings.pseudopotentials, 
                    input_data=settings.espresso_settings_dict,
                    filename=filename,
                    kpts= settings.kpoints[1] if 'gamma' not in settings.kpoints[0] else None, koffset=settings.kpoints[2] if 'gamma' not in settings.kpoints[0] else None)
        if(settings.fixed_layers_slab): fixed_slab = slab.get_atoms_by_layers(settings.fixed_layers_slab)
        else: fixed_slab = settings.fixed_indices_slab
        calc.set_fixed_atoms(fixed_slab, slab.reindex_map, settings.fixed_indices_mol, mol.reindex_map, slab.natoms, mol.natoms, settings.fix_slab_xyz, settings.fix_mol_xyz)
        calc.set_system_flags(settings.starting_mag, settings.flags_i)
        calc.write_input(all_mol_on_slab_configs_ase[i])

    #print(pwi_names)
    launch_jobs(jobscript=settings.jobscript, pwi_list=pwi_names, outdirs=relax_outdir, jobname_title='rel')


def saveas(which : str, i_or_f : str, saveas_format : str):

    if i_or_f == 'i':
        pw = 'pwi'
    elif i_or_f == 'f':
        pw = 'pwo'
    else: 
        raise RuntimeError("Wrong arguments: passed '{0} {1}', expected 'screening i/f' or 'relax i/f'".format(which, i_or_f))
    if which != 'screening' and which != 'relax':
        raise RuntimeError("Wrong argument: passed '{0}', expected 'screening' or 'relax'".format(which))

    print('Reading files...')
    pw_list=glob.glob(pw_files_prefix+which+'_*.'+pw)
    configs = [(read(file) if pw == 'pwi' else read(file, results_required=False)) for file in pw_list]
    print('All files read.')

    print("Saving {0} files to {1} format...".format(which, saveas_format))
    folder = saveas_format+'/'+which+'/'
    if not os.path.exists(saveas_format):
        os.mkdir(saveas_format)
    if not os.path.exists(folder):
        os.mkdir(folder)    

    for i, config in enumerate(configs):
        write(folder+pw_list[i].replace(pw, saveas_format), config)

    print("Files saved to {0}".format(saveas_format+'/'+which) )
    

