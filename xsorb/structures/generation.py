# from multipledispatch import dispatch

# class Shape:
#     @dispatch(int)
#     def __init__(self, length):
#         print("int")

#     @dispatch(float)
#     def __init__(self, length):
#         print("float")

#     @dispatch()
#     def __init__(self):
#         print("No argument")


# if __name__ == "__main__":
#     line = Shape(3)
#     curve = Shape(5.5)
#     point = Shape()



def generate_adsorption_structures(self, molecule : Atoms, 
                                    adsites : list, 
                                    z_distance_from_site : float,
                                    min_z_distance_from_surf : float,
                                    adsites_labels : list,
                                    rotation_label : str,
                                    connected_adsites : dict = None,
                                    surrounding_exclude_main : bool = False,
                                    mol_before_slab : bool = False):
    '''
    Returns the adsorption structures obtained by placing the molecule on all the adsites,
    ensuring that the molecule is not too close to the surface

    Args:
    - molecule: ase Atoms for the molecule to be placed in all the adsorption sites
    - adsites: list of cartesian coordinates of the adsorption sites
    - z_distance_from_site: target distance between reference atom and site.
    - min_z_distance_from_surf: minimum vertical distance between ANY atom of the molecule and ANY atom of the surface
    - adsites_labels: list of labels of the adsorption sites, used to create the final full labels of the adsorption configurations
    - rotation_label: label of the rotation of the molecule, used to create the final full labels of the adsorption configurations
    - mol_before_slab: decide wether to put the molecule before or after the slab in the list of atoms
    '''

    adsorption_structures = []
    full_labels = []

    if not connected_adsites:
        for coords, site_label in zip(adsites, adsites_labels):
            mol = molecule.copy()

            #place the molecule in the adsorption site, then translate it upwards by the target height
            mol.translate(coords) 
            mol.translate([0,0,z_distance_from_site])
            final_deltaz = z_distance_from_site

            #Check for min_z_distance_from_surf and translate accordingly
            dz = mindistance_deltaz(self.slab_ase, mol, min_z_distance_from_surf)
            if dz:
                mol.translate( [0, 0, dz] )
                #print("The molecule was translated further by {0} to enforce minimum distance.".format(dz))
                final_deltaz += dz
                
            struct = mol + self.slab_ase if mol_before_slab else self.slab_ase + mol
            struct.cell = self.slab_ase.cell
            adsorption_structures.append(struct)
            full_labels.append(rotation_label+site_label+'{:.3f}'.format(final_deltaz))
    
    else:
        yrot = int(float(rotation_label.split(',')[1]))
        for main_site_num, main_site_coords, main_site_label in zip(range(len(adsites)), adsites, adsites_labels):
            
            # create the rotations for when the molecule is horizontal,
            #orienting it towards the nearest neighbors
            
            if yrot != 90:
                for rel_ads in connected_adsites[main_site_num]:

                    if surrounding_exclude_main and rel_ads['duplicate_main']:
                        continue

                    mol = molecule.copy()

                    #rotate the molecule towards the surrounding site
                    main_to_rel_axis = (rel_ads['position'] - main_site_coords)
                    main_to_rel_axis[2] = 0 #project on xy plane
                    mol.rotate('x', main_to_rel_axis)

                    #place the molecule in the adsorption site, then translate it upwards by the target height
                    mol.translate(main_site_coords) 
                    mol.translate([0,0,z_distance_from_site])
                    final_deltaz = z_distance_from_site

                    #Check for min_z_distance_from_surf and translate accordingly
                    dz = mindistance_deltaz(self.slab_ase, mol, min_z_distance_from_surf)
                    if dz:
                        mol.translate( [0, 0, dz] )
                        #print("The molecule was translated further by {0} to enforce minimum distance.".format(dz))
                        final_deltaz += dz
                        
                    struct = mol + self.slab_ase if mol_before_slab else self.slab_ase + mol
                    struct.cell = self.slab_ase.cell
                    adsorption_structures.append(struct)
                    rotation_label = rotation_label.split(',')
                    rotation_label[2] = 'to_'+rel_ads['label'].split()[0]
                    rotation_label = ','.join(rotation_label)
                    full_labels.append(rotation_label+main_site_label+'{:.3f}'.format(final_deltaz))
            
            else: # vertical molecule
                rel_ads_positions = [rel_ads['position'] for rel_ads in connected_adsites[main_site_num] \
                                        if not rel_ads['duplicate_surrounding'] and not rel_ads['duplicate_main']]
                rel_ads_labels = [rel_ads['label'] for rel_ads in connected_adsites[main_site_num] \
                                    if not rel_ads['duplicate_surrounding'] and not rel_ads['duplicate_main']]
                
                for coords, site_label in zip([adsites[main_site_num]]+rel_ads_positions, [adsites_labels[main_site_num]]+rel_ads_labels):
                    mol = molecule.copy()

                    #place the molecule in the adsorption site, then translate it upwards by the target height
                    mol.translate(coords) 
                    mol.translate([0,0,z_distance_from_site])
                    final_deltaz = z_distance_from_site

                    #Check for min_z_distance_from_surf and translate accordingly
                    dz = mindistance_deltaz(self.slab_ase, mol, min_z_distance_from_surf)
                    if dz:
                        mol.translate( [0, 0, dz] )
                        #print("The molecule was translated further by {0} to enforce minimum distance.".format(dz))
                        final_deltaz += dz
                        
                    struct = mol + self.slab_ase if mol_before_slab else self.slab_ase + mol
                    struct.cell = self.slab_ase.cell
                    adsorption_structures.append(struct)
                    full_labels.append(rotation_label+site_label+'{:.3f}'.format(final_deltaz))



    return adsorption_structures, full_labels


def closest_pair(slab : Atoms, mol: Atoms):
    '''
    Returns the indices of the closest pair of (slab, molecule) atoms, and their distance

    Args:
    - slab: Atoms object for the slab
    - mol: Atoms object for the molecule 
    '''

    dist_matrix = get_distances(mol.positions, slab.positions, slab.cell, pbc=True)[1]

    i_mol, j_slab = np.unravel_index(np.argmin(dist_matrix), dist_matrix.shape)
    mindist = dist_matrix[i_mol, j_slab]

    return i_mol, j_slab, mindist


def mindistance_deltaz(slab : Atoms, mol: Atoms, min_z_distance_from_surf : float):
    '''
    Returns the vertical translation required so that the closest mol-slab atom pair is at least
    min_z_distance_from_surf apart (or half the sum of the covalent radii) along z

    Args:
    - slab: Atoms object for the slab
    - mol: Atoms object for the molecule
    - min_distance: minimum required distance along z
    '''

    dz_tot = 0
    molcopy = mol.copy()
    while(True):
        #First, find the closest slab-mol atoms pair
        i_mol, j_slab, _ = closest_pair(slab, molcopy)

        #Find the z coordinates of the closest atoms pair, and half their covalent distance
        half_covalent_distance = 0.5 * (covalent_radii[atomic_numbers[slab[j_slab].symbol]] \
                                    + covalent_radii[atomic_numbers[molcopy[i_mol].symbol]])
        zmol = molcopy[i_mol].position[2]
        zslab = slab[j_slab].position[2]

        #Calculate the distance required to enforce the minimum distance
        necessary_min_z_dist = max(min_z_distance_from_surf, half_covalent_distance)
        if(zmol < zslab + necessary_min_z_dist):
            dz = zslab + necessary_min_z_dist - zmol
            dz_tot += dz
            molcopy.translate([0,0,dz])
        else:
            break
    
    return dz_tot


def slab_mol_bonds(slab : Atoms, mol: Atoms):
    '''
    Returns the list of the bonds between the molecule and the slab.

    Args:
    - slab: Atoms object for the slab
    - mol: Atoms object for the molecule
    '''

    atoms = slab+mol
    cutoffs = natural_cutoffs(atoms, mult=1.15)
    nl = NeighborList(cutoffs, skin=0, self_interaction=False, bothways=True)
    nl.update(atoms)
    cm = nl.get_connectivity_matrix()

    bonds_list = []
    for i in range(len(slab)):
        for j in range(len(mol)):
            if cm[i, len(slab)+j]:
                bonds_list.append('{0}{1}-{2}{3}'.format(mol.symbols[j], j, slab.symbols[i], i))
    
    return bonds_list
