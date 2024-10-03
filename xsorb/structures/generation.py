from dataclasses import dataclass, asdict

from ase import Atoms
from ase.data import atomic_numbers, covalent_radii

from xsorb.structures.utils import closest_pair
from xsorb.settings import Settings
from xsorb.structures.molecule import Molecule
from xsorb.structures.slab import Slab, AdsorptionSite
from xsorb.structures.molecule import MoleculeRotation


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


@dataclass
class AdsorptionStructure:
    atoms: Atoms
    adsite: AdsorptionSite
    mol_rot: MoleculeRotation
    dz : float

    def to_dataframe_row(self):
        return {
            # 'adsite': self.adsite.label,
            # 'mol_rot': self.mol_rot.to_dict(),
            # 'dz': self.dz
        }


class AdsorptionStructuresGenerator:

    def __init__(self, settings: Settings, VERBOSE : bool = False) -> None:
    
        self.settings = settings
        self.VERBOSE = VERBOSE

        if self.settings.structure.adsorption_sites.mode == 'coord_number'\
            and self.settings.structure.adsorption_sites.coord_number_params.include_surrounding_sites:
            self.rot_mode = 'surrounding'
        else:
            self.rot_mode = 'standard'

        self.mol_before_slab = self.settings.structure.misc.mol_before_slab
    
        #Slab import from file
        if self.VERBOSE: print('Loading slab...')
        self.slab = Slab(slab_filename=settings.input.slab_filename,
                    surface_sites_height=settings.structure.adsorption_sites.surface_height, 
                    layers_threshold=settings.structure.constraints.layers_height, 
                    fixed_layers_slab=settings.structure.constraints.fixed_layers_slab, 
                    fixed_indices_slab=settings.structure.constraints.fixed_indices_slab, 
                    fix_slab_xyz=settings.structure.constraints.fix_slab_xyz,
                    sort_atoms_by_z=settings.structure.misc.sort_atoms_by_z,
                    translate_slab_from_below_cell_bottom=settings.structure.misc.translate_slab)
        if self.VERBOSE: print('Slab loaded.')


        #Molecule import from file
        if self.VERBOSE: print('Loading molecule...')
        self.mol = Molecule(molecule_filename=settings.input.molecule_filename,
                    atom_index=settings.structure.molecule.selected_atom_index,
                    molecule_axis_atoms=settings.structure.molecule.molecule_axis["values"] \
                        if settings.structure.molecule.molecule_axis["mode"] == "atom_indices" else None, 
                    axis_vector=settings.structure.molecule.molecule_axis["values"] \
                        if settings.structure.molecule.molecule_axis["mode"] == "vector" else None,  
                    fixed_indices_mol=settings.structure.constraints.fixed_indices_mol, 
                    fix_mol_xyz=settings.structure.constraints.fix_mol_xyz)
        if self.VERBOSE: print('Molecule loaded.')



    def generate_molecule_rotations(self,
                                     rot_mode : str = 'standard',
                                     adsite : AdsorptionSite | None = None, 
                                     SAVEFIG : bool = False):

        if rot_mode == 'standard':
            if hasattr(self, 'molecule_rotations'):
                return self.molecule_rotations
            z_rot_angles = self.settings.structure.molecule.z_rot_angles
            surrounding_exclude_main = False
            VERBOSE = self.VERBOSE
        else:
            if adsite is None:
                raise ValueError("adsite must be provided for surrounding mode")
            z_rot_angles = adsite.surrounding_sites
            surrounding_exclude_main = \
                self.settings.structure.adsorption_sites.coord_number_params.surrounding_exclude_main
            VERBOSE = False

        molecule_rotations = self.mol.generate_molecule_rotations(
            x_rot_angles=self.settings.structure.molecule.x_rot_angles, 
            y_rot_angles=self.settings.structure.molecule.y_rot_angles,
            z_rot_angles=z_rot_angles, 
            vert_angles_list=self.settings.structure.molecule.vertical_angles,
            individual_rotations=self.settings.structure.molecule.individual_rotations,
            surrounding_exclude_main=surrounding_exclude_main,
            SAVE_IMAGE=SAVEFIG,
            VERBOSE=VERBOSE)

        if rot_mode == 'standard': 
            #in case of standard, store the rotations to avoid calling the function multiple times
            self.molecule_rotations = molecule_rotations

        return molecule_rotations


    def put_together_slab_and_mol(self, adsite : AdsorptionSite, 
                                  mol_rot : MoleculeRotation):
            mol = mol_rot.atoms.copy()
            
            #place the molecule in the adsorption site, then translate it upwards by the target height
            mol.translate(adsite.coords)

            dz = self.settings.structure.molecule.target_distance
            mol.translate([0,0,dz])

            #Check for min_z_distance_from_surf and translate accordingly
            further_transl = mindistance_deltaz(self.slab.slab_ase, 
                                                mol, 
                                                self.settings.structure.molecule.min_distance)
            if further_transl:
                mol.translate( [0, 0, further_transl] )
                #print("The molecule was translated further by {0} to enforce minimum distance.".format(dz))
                dz += further_transl

            struct : Atoms = mol + self.slab.slab_ase if self.mol_before_slab else self.slab.slab_ase + mol
            struct.cell = self.slab.slab_ase.cell

            return AdsorptionStructure(atoms=struct, adsite=adsite, mol_rot=mol_rot, dz=dz)


    def get_structures_for_vertical_surrounding_sites(self, adsite : AdsorptionSite):
        
        #put the molecule vertically on the surrounding sites

        adsorption_structures_site : list[AdsorptionStructure] = []

        for surr_site in adsite.surrounding_sites:

            if surr_site.duplicate_main or surr_site.duplicate_surrounding:
                continue
            
            #promote surr_site to AdsorptionSite
            surr_site = AdsorptionSite(label=surr_site.label, 
                                        coords=surr_site.coords, 
                                        type='ontop',
                                        info='',
                                        unique_id=surr_site.unique_id)
            
            molecule_rotations = self.generate_molecule_rotations(rot_mode='standard')
            for mol_rot in molecule_rotations:
                #if y_angle != +/- 90, continue
                if mol_rot.yrot not in ['90', '-90']:
                    continue
                structure = self.put_together_slab_and_mol(adsite=surr_site,mol_rot=mol_rot)
                adsorption_structures_site.append(structure)


        return adsorption_structures_site


    def generate_adsorption_structures(self, SAVEFIG : bool = False):
        '''
        Generate all adsorption structures considering the combinations of 
        molecular rotations and adsorption sites.

        Returns a list of  ASE atoms with the configurations, 
        and a list of info (csv format) of each configuration.

        Args:
        - settings: Settings object, containing the filenames of molecule and slab, and all the parameters for the rotations of the molecule and the indentifications of adsorption sites
        - SAVEFIG: save an image of the adsorption sites and of the molecular rotations.
        ''' 


        #Find adsorption sites and labels (site type and x,y coords.)
        mode_params = {'high_symmetry': self.settings.structure.adsorption_sites.high_symmetry_params, 
                 'coord_number': self.settings.structure.adsorption_sites.coord_number_params,}
        mode = self.settings.structure.adsorption_sites.mode    
        if mode not in mode_params:
            raise ValueError(f"mode must be one of {mode_params.keys()}")


        adsites = self.slab.find_adsorption_sites(
            mode=mode,
            **asdict(mode_params[mode]),
            selected_sites=self.settings.structure.adsorption_sites.selected_sites,
            SAVE_IMAGE=SAVEFIG,
            VERBOSE=self.VERBOSE)


        #Adsorption of molecule on all adsorption sites for all molecule orientations
        if self.VERBOSE: print('Generating adsorption structures...') 

        adsorption_structures : list[AdsorptionStructure] = []
        for adsite in adsites:
            molecule_rotations = self.generate_molecule_rotations(rot_mode=self.rot_mode, adsite=adsite,SAVEFIG=SAVEFIG)
            for mol_rot in molecule_rotations:
                structure = self.put_together_slab_and_mol(adsite=adsite,mol_rot=mol_rot)
                adsorption_structures.append(structure)


            #handle the case of vertical molecule on surrounding sites
            if adsite.surrounding_sites is not None:
                adsorption_structures.extend(self.get_structures_for_vertical_surrounding_sites(adsite))


        if self.VERBOSE: print('Adsorption structures generated.')

        return adsorption_structures