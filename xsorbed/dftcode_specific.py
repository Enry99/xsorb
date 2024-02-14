"""
@author: Enrico Pedretti

DFT CODE-specific functions (all others must be code-agnostic)

"""

import shutil, glob, os, json
from ase.units import create_units
from ase.constraints import FixScaled
from ase.calculators.espresso import Espresso
from ase.calculators.vasp import Vasp
from xsorbed.common_definitions import *

#TODO: check that when reading magmoms from input, the order is then changed correctly after resorting the poscar during write_inputs
#TODO: Check that ase-sort.dat is used to sort the files back in the initial order when reading the outputs.
#It does not work when reading poscar, it should work with outcar/vasprun.xml

SUPPORTED_PROGRAMS = ['VASP', 'ESPRESSO']

HYBRID_SCREENING_THRESHOLDS = {
    'VASP' : [-0.5], # ~ -2e-2 Ry/Bohr
    'ESPRESSO' : [5e-3, 5e-2]
}


UNITS_TO_EV_FACTOR = {
    'VASP' : 1,
    'ESPRESSO': create_units('2006')['Rydberg']
}

COMPLETION_STRINGS = {
    'RELAX_COMPLETED': {
        #'VASP' : 'reached required accuracy - stopping structural energy minimisation', #in OUTCAR
        'VASP' : 'finalpos', #in vasprun.xml
        'ESPRESSO': 'Final energy'
    },

    'SCF_NONCONVERGED': {
        'VASP': 'abcdefgxyz', #TODO: vasp does not stop the relax if one scf does not converge, so not necessary
        'ESPRESSO': 'convergence NOT achieved'
    }
}

OUT_FILE_PATHS = {
    'SCREENING': {
        'VASP': screening_outdir+'/{0}/vasprun.xml',
        'ESPRESSO': 'screening_{0}.pwo',
    },

    'RELAX': {
        'VASP': relax_outdir+'/{0}/vasprun.xml',
        'ESPRESSO': 'relax_{0}.pwo',
    }
}

IN_FILE_PATHS = {
    'SCREENING': {
        'VASP': screening_outdir+'/{0}/POSCAR',
        'ESPRESSO': 'screening_{0}.pwi',
    },

    'RELAX': {
        'VASP': relax_outdir+'/{0}/POSCAR',
        'ESPRESSO': 'relax_{0}.pwi',
    }
}

SBATCH_POSTFIX = {
    'SCREENING': {
        'VASP': '',
        'ESPRESSO': '{0}/screening_{1}.pwi {0}/screening_{1}.pwo',
    },

    'RELAX': {
        'VASP': '',
        'ESPRESSO': '{0}/screening_{1}.pwi {0}/screening_{1}.pwo',
    }
}


FRAGMENTS_IN_FILE_PATHS = {
    'VASP': 'fragments/{0}/POSCAR',
    'ESPRESSO': 'fragments/{0}/{0}.pwi'
}

FRAGMENTS_OUT_FILE_PATHS = {
    'VASP': 'fragments/{0}/vasprun.xml',
    'ESPRESSO': 'fragments/{0}/{0}.pwo'
}

SBATCH_POSTFIX_FRAGS = {
    'VASP': '',
    'ESPRESSO': '{0}.pwi {0}.pwo',
}


def parse_espresso_settings(block_str_list : list):
    
    #NOTE 1: The blocks CELL_PARAMETERS ATOMIC_POSITIONS ATOMIC_SPECIES must NOT be included in input file, 
    #as they are read from the input structures
    #NOTE 2: This code does not yet support the following Espresso blocks:
    #OCCUPATIONS, CONSTRAINTS, ATOMIC_VELOCITIES, ATOMIC_FORCES, ADDITIONAL_K_POINTS, SOLVENTS, HUBBARD
    
    from ase.io.espresso import read_fortran_namelist    
    
    # parse namelist section and extract remaining lines
    dftprogram_settings_dict, card_lines = read_fortran_namelist(block_str_list)

    #parse ATOMIC_SPECIES and K_POINTS
    for i, line in enumerate(card_lines):
        if('ATOMIC_SPECIES' in line.upper()):
            atomic_species_index = i
        if('K_POINTS' in line.upper()):
            k_points_index = i

    #ATOMIC_SPECIES
    i = atomic_species_index+1

    dftprogram_settings_dict['pseudopotentials'] = {}
    while i < len(card_lines):
        line = card_lines[i]
        
        if ('ATOMIC_POSITIONS' in line.upper() or 
        'K_POINTS' in line.upper() or 
        'ADDITIONAL_K_POINTS' in line.upper() or 
        'CELL_PARAMETERS' in line.upper() or 
        'CONSTRAINTS' in line.upper() or 
        'OCCUPATIONS' in line.upper() or 
        'ATOMIC_VELOCITIES' in line.upper() or 
        'ATOMIC_FORCES' in line.upper() or 
        'SOLVENTS' in line.upper() or 
        'HUBBARD' in line.upper()): break       
        
        element, mass, pseudo = line.split()
        dftprogram_settings_dict['pseudopotentials'].update({element : pseudo})
        i+=1

    #K_POINTS
    if 'gamma' in card_lines[k_points_index].split()[1].strip().lower():
        dftprogram_settings_dict['kpts'] = None
        dftprogram_settings_dict['koffset'] = None
    else:
        line = card_lines[k_points_index+1]
        dftprogram_settings_dict['kpts'] = line.split()[:3]
        dftprogram_settings_dict['koffset'] = line.split()[3:]

    return dftprogram_settings_dict


def parse_vasp_settings(block_str_list : list):

    dftprogram_settings_dict = {}
    cards = {'GENERAL' : [], 'INCAR' : [], 'KPOINTS': []}
    
    in_card = False
    for line in block_str_list:
        
        #begin card
        if line.strip()[0] == '&':
            block_name = line.split('&')[1].strip()
            in_card = True
            continue

        #end card
        if line.strip()=='/':
            in_card = False

        #within card
        if in_card: cards[block_name.upper()].append(line)  
    

    for line in cards['GENERAL']:
        (key, val) = line.split('=')
        key = key.strip()
        val = val.strip()
        val = val.strip("'")
        val = val.strip('"')   
        key = key.lower()     
        
        if key == 'vasp_pseudo_setups':
            val = json.loads(val)
        dftprogram_settings_dict[key] = val
        
    if cards['INCAR']: dftprogram_settings_dict["incar_string"] = '\n'.join(cards['INCAR']) #use as an untouched single string
    if cards['KPOINTS']: dftprogram_settings_dict["kpoints_string"] = '\n'.join(cards['KPOINTS']) #use as an untouched single string

    if 'pymatgen_set' not in dftprogram_settings_dict and 'incar_string' not in dftprogram_settings_dict:
        pass
        #raise RuntimeError("")
    
    if "vasp_xc_functional" not in dftprogram_settings_dict:
        dftprogram_settings_dict["vasp_xc_functional"] = 'PBE'
    
    if "vasp_pseudo_setups" not in dftprogram_settings_dict:
        dftprogram_settings_dict["vasp_pseudo_setups"] = {"base": "recommended"} 


    return dftprogram_settings_dict
    

def get_dftprogram_settings(program : str, block_str_list : list):
    
    if program == 'ESPRESSO':
        return parse_espresso_settings(block_str_list)
    elif program == 'VASP':
        return parse_vasp_settings(block_str_list)
    else:
        raise RuntimeError('Program not recognized')




def override_settings(settings, calc_type : str):
    '''
    Overrides the DFT code settings with specific options for the preliminary screening or the final relaxation

    Args:
    - settings: Settings object, containing the dft code parameters
    - calc_type: 'SCREENING' or 'RELAX'
    '''
    
    if calc_type == 'SCREENING':
        if settings.program == 'ESPRESSO':

            settings.dftprogram_settings_dict['control'].update({'calculation' : 'relax' })
            settings.dftprogram_settings_dict['control'].update({'restart_mode' : 'from_scratch'})
            settings.dftprogram_settings_dict['control'].update({'etot_conv_thr' : settings.screening_conv_thr[0]})
            settings.dftprogram_settings_dict['control'].update({'forc_conv_thr' : settings.screening_conv_thr[1]})
            settings.dftprogram_settings_dict['control'].update({'outdir' : 'WORK'}) #TODO: lasciarlo scegliere all'utente

            if 'ions' not in settings.dftprogram_settings_dict: settings.dftprogram_settings_dict['ions'] = {}
            settings.dftprogram_settings_dict['ions'].update({'upscale': 1})


        elif settings.program == 'VASP':

            #fix if user forgot to put the correct IBRION for relax, and add EDIFFG for screening.
            #EDIFFG is put here so that in any case this will be the final value, even if the user had 
            #specified a different value explicitly in the &INCAR instead of using one of the RelaxSets
            if 'incar_string' in settings.dftprogram_settings_dict:
                s = settings.dftprogram_settings_dict['incar_string'].split('\n')
            else: s = []

            missing_ibrion = True
            missing_ediffg = True
            for i, line in enumerate(s):
                if 'EDIFFG' in line:
                    s[i] = f'EDIFFG = {settings.screening_conv_thr[0]}'
                    missing_ediffg = False
                if 'IBRION' in line: 
                    missing_ibrion = False
            if missing_ediffg: s.append(f'EDIFFG = {settings.screening_conv_thr[0]}')
            if missing_ibrion: s.append('IBRION = 2')

            settings.dftprogram_settings_dict['incar_string'] = '\n'.join(s)


    elif calc_type == 'RELAX':
        if settings.program == 'ESPRESSO':
            settings.dftprogram_settings_dict['control'].update({'calculation' : 'relax'})
            settings.dftprogram_settings_dict['control'].update({'restart_mode' : 'from_scratch'})
            settings.dftprogram_settings_dict['control'].update({'outdir' : 'WORK'}) #TODO: lasciarlo scegliere all'utente
            if 'ions' not in settings.dftprogram_settings_dict: settings.dftprogram_settings_dict['ions'] = {}
            #settings.dftprogram_settings_dict['ions'].update({'ion_dynamics': settings.ion_dynamics})
        
        elif settings.program == 'VASP':
            
            #fix if user forgot to put the correct IBRION for relax
            if 'incar_string' in settings.dftprogram_settings_dict and "pymatgen_set" not in settings.dftprogram_settings_dict:
                s = settings.dftprogram_settings_dict['incar_string'].split('\n')
            
                missing_ibrion = True
                for i, line in enumerate(s):
                    if 'IBRION' in line: missing_ibrion = False
                if missing_ibrion: s.append('IBRION = 2')

                settings.dftprogram_settings_dict['incar_string'] = '\n'.join(s)

def override_settings_isolated_fragment(settings, natoms_mol : int, manual_dft_override : dict = None):
    '''
    Overrides the DFT code settings with specific options for the isolated fragment calculation

    Args:
    - settings: Settings object, containing the dft code parameters
    - natoms_mol: number of atoms in the molecule
    - manual_dft_override: additional settings fragment-specific, read from fragments.json
    '''
    if settings.program == 'ESPRESSO':

        settings.dftprogram_settings_dict['control'].update({'outdir' : 'WORK'}) #TODO: lasciarlo scegliere all'utente
        settings.dftprogram_settings_dict['system'].update({'nosym' : True})
        settings.dftprogram_settings_dict['system'].update({'starting_magnetization(1)' : 1.0})        
        settings.dftprogram_settings_dict['electrons'].update({'mixing_beta' : 0.1})

        if manual_dft_override is not None:
            if 'SYSTEM' in manual_dft_override:
                for k, v in manual_dft_override['SYSTEM'].items():
                    settings.dftprogram_settings_dict['system'][k] = v
            if 'ELECTRONS' in manual_dft_override:
                for k, v in manual_dft_override['ELECTRONS'].items():
                    settings.dftprogram_settings_dict['electrons'][k] = v


    elif settings.program == 'VASP':
        if 'incar_string' in settings.dftprogram_settings_dict:
            s = settings.dftprogram_settings_dict['incar_string'].split('\n')
        else: s = []

        missing_isym = True
        missing_magmom = True
        for i, line in enumerate(s):
            
            if 'ISYM' in line:
                if manual_dft_override and 'ISYM' in manual_dft_override:
                    s[i] = f'ISYM = {manual_dft_override["ISYM"].strip()}'
                    manual_dft_override.pop('ISYM')
                else:
                    s[i] = 'ISYM = 0'          
                missing_isym = False

            if 'MAGMOM' in line:
                if manual_dft_override and 'MAGMOM' in manual_dft_override:
                    s[i] = f'MAGMOM = {manual_dft_override["MAGMOM"].strip()}'
                    manual_dft_override.pop('MAGMOM')
                else:
                    s[i] = f'MAGMOM = {natoms_mol}*0.6'
                missing_magmom = False
            
            key = line.split()[0].strip()
            if key in manual_dft_override:
                s[i] = f'{key} = {manual_dft_override[key].strip()}'
                manual_dft_override.pop(key)
                
        if missing_isym: s.append('ISYM = 0')
        if missing_magmom: s.append(f'MAGMOM = {natoms_mol}*0.6')

        if manual_dft_override is not None:
            for key in manual_dft_override:
                s.append(f'{key} = {manual_dft_override[key].strip()}')


        settings.dftprogram_settings_dict['incar_string'] = '\n'.join(s)     

def override_settings_adsorbed_fragment(program : str, dft_section : list, dft_settings_dict : dict = None):
    '''
    Overrides the DFT code settings with specific options for the adsorbed fragments calculations
    '''

    dft_newlines = []
    
    if program == 'ESPRESSO':
    
        card = None
        for line in dft_section:
            
            if '&' in line: 
                card = line.split('&')[1].strip()
                dft_newlines.append(line)
                continue
            
            if card and card in dft_settings_dict:
                found = False
                for key, val in dft_settings_dict[card].items():
                    if key in line: 
                        dft_newlines.append(f'   {key} = {val}\n')
                        dft_settings_dict[card].pop(key)
                        found = True
                        break
                if found: continue
                
            if card and line.strip()=='/': 
                for key, val in dft_settings_dict[card].items():
                    dft_newlines.append(f'   {key} = {val}\n')
                dft_newlines.append(line)
                card = None


    elif program == 'VASP':
        
        in_incar = False
        for line in dft_section:
            
            if '&INCAR' in line: 
                in_incar = True
                dft_newlines.append(line)
                continue
            
            if in_incar:
                found = False
                for key, val in dft_settings_dict.items():
                    if key in line: 
                        dft_newlines.append(f'   {key} = {val}\n')
                        dft_settings_dict.pop(key)
                        found = True
                        break
                if found: continue
                
            if in_incar and line.strip()=='/': 
                for key, val in dft_settings_dict.items():
                    dft_newlines.append(f'   {key} = {val}\n')
                dft_newlines.append(line)
                card = None
        

    return dft_newlines    



def Calculator(settings, label : str, atoms, directory : str):
    
    if settings.program == 'ESPRESSO':
        return Espresso(label = label,
                        pseudopotentials=settings.dftprogram_settings_dict['pseudopotentials'], 
                        input_data=settings.dftprogram_settings_dict)
                            
    elif settings.program == 'VASP':
        
        if "pymatgen_set" in settings.dftprogram_settings_dict:
            from pymatgen.io.vasp.sets import MPRelaxSet, MPMetalRelaxSet, MPScanRelaxSet, MPHSERelaxSet, MITRelaxSet
            from pymatgen.io.ase import AseAtomsAdaptor

            if settings.dftprogram_settings_dict["pymatgen_set"].lower() == 'MPRelaxSet'.lower(): RelaxSet = MPRelaxSet
            elif settings.dftprogram_settings_dict["pymatgen_set"].lower() == 'MPMetalRelaxSet'.lower(): RelaxSet = MPMetalRelaxSet
            elif settings.dftprogram_settings_dict["pymatgen_set"].lower() == 'MPScanRelaxSet'.lower(): RelaxSet = MPScanRelaxSet
            elif settings.dftprogram_settings_dict["pymatgen_set"].lower() == 'MPHSERelaxSet'.lower(): RelaxSet = MPHSERelaxSet            
            elif settings.pymatgen_set.lower() == 'MITRelaxSet'.lower(): RelaxSet = MITRelaxSet
            else: raise ValueError('Pymatgen preset not recognized.')

            atomscopy = atoms.copy()
            del atomscopy.constraints #to suppress the warning about constraints not supported in pymatgen
            relax = RelaxSet(AseAtomsAdaptor.get_structure(atomscopy))
            preset_incar_settings = {k.lower(): v for k, v in relax.incar.as_dict().items()} 
            if 'isif' in preset_incar_settings: preset_incar_settings['isif'] = 2 #OVERRIDE: we do not want a vc-relax, unless explicitly specified in &INCAR
            preset_incar_settings.pop('@module')
            preset_incar_settings.pop('@class')
            relax.kpoints.write_file('_temp_kpts_')

        preset_incar_settings["xc"] = settings.dftprogram_settings_dict["vasp_xc_functional"]

        adjust_constraints(atoms, 'VASP')
        os.environ["VASP_PP_PATH"] = settings.dftprogram_settings_dict["vasp_pp_path"]
        calc = Vasp(directory=directory, 
                    setups=settings.dftprogram_settings_dict["vasp_pseudo_setups"], 
                    **preset_incar_settings) #set here the default values from pymat recommended


        #write user-defined settings to string, to be parsed by ASE
        if "incar_string" in settings.dftprogram_settings_dict:
            with open('_temp_incar_', 'w') as f:
                f.write(settings.dftprogram_settings_dict["incar_string"])
            calc.read_incar('_temp_incar_')
            os.remove('_temp_incar_') 

        if "kpoints_string" in settings.dftprogram_settings_dict:
            with open('_temp_kpts_', 'w') as f:
                f.write(settings.dftprogram_settings_dict["kpoints_string"])

        if "kpoints_string" in settings.dftprogram_settings_dict or "pymatgen_set" in settings.dftprogram_settings_dict:
            calc.read_kpoints('_temp_kpts_') 
            os.remove('_temp_kpts_')
        
        return calc



def adjust_constraints(atoms, program : str):
    
    if program == 'VASP':
        c = [FixScaled(atoms.cell, constr.a, ~constr.mask) for constr in atoms.constraints] #atoms.constraints are FixCartesian
        atoms.set_constraint(c)


def edit_files_for_restart(program : str, calc_type : str, indices : list):
    '''
    Edit the input files, setting the correct flags for restart.

    Args:
    - program: DFT program. Possible values: 'ESPRESSO' or 'VASP'
    - calc_type: 'SCREENING' or 'RELAX'
    - indices_list: indices of the calculations
    '''

    for index in indices:
        if program == 'ESPRESSO':
            with open(IN_FILE_PATHS[calc_type][program].format(index), 'r') as f:
                lines = f.readlines()
                for i, line in enumerate(lines):
                    if 'from_scratch' in line:
                        lines[i] = lines[i].replace('from_scratch','restart')
                        break
                f.writelines(lines)
            with open(IN_FILE_PATHS[calc_type][program].format(index), 'w') as f:
                f.writelines(lines)

        elif program == 'VASP':
            poscar = IN_FILE_PATHS[calc_type][program].format(index)
            contcar = poscar.replace('POSCAR', 'CONTCAR')
            incar = poscar.replace('POSCAR', 'INCAR')
            outcar = poscar.replace('POSCAR', 'OUTCAR')
            vasprun = poscar.replace('POSCAR', 'vasprun.xml')
            oszicar = poscar.replace('POSCAR', 'OSZICAR')
            
            with open(incar, 'r') as f:
                lines = f.readlines()
                istart_found = False
                for i, line in enumerate(lines):
                    if 'ISTART' in line:
                        lines[i] = 'ISTART = 1\n'
                        istart_found = True
                        break
                if not istart_found: lines.append('ISTART = 1\n')
                f.writelines(lines)

            with open(incar, 'w') as f:
                f.writelines(lines)

            #copy files for the first part of the relaxation, to avoid overwriting
            last_i = len( glob.glob( outcar.replace('OUTCAR', 'OUTCAR*') ) ) - 1
            shutil.copyfile(outcar, outcar.replace('OUTCAR', f'OUTCAR_{last_i}'))
            shutil.copyfile(vasprun, vasprun.replace('vasprun.xml', f'vasprun_{last_i}.xml'))
            shutil.copyfile(poscar, poscar.replace('POSCAR', f'POSCAR_{last_i}'))
            shutil.copyfile(oszicar, oszicar.replace('OSZICAR', f'OSZICAR_{last_i}'))
           
            #copy contcar to poscar to restart
            shutil.copyfile(contcar, poscar)
