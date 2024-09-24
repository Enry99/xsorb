from dataclasses import dataclass, field
from typing import Optional
from dacite import from_dict, Config
from xsorb.dft_codes.definitions import HYBRID_SCREENING_THRESHOLDS



@dataclass
class EspressoSettings:
    pwi_path: str
    etot_conv_thr_screening: float = HYBRID_SCREENING_THRESHOLDS['ESPRESSO'][0]
    forc_conv_thr_screening: float = HYBRID_SCREENING_THRESHOLDS['ESPRESSO'][1]


def read_espresso_settings(input_settings_dict : dict):
    
    #NOTE 1: The blocks CELL_PARAMETERS ATOMIC_POSITIONS ATOMIC_SPECIES must NOT be included in input file, 
    #as they are read from the input structures
    #NOTE 2: This code does not yet support the following Espresso blocks:
    #OCCUPATIONS, CONSTRAINTS, ATOMIC_VELOCITIES, ATOMIC_FORCES, ADDITIONAL_K_POINTS, SOLVENTS
    

    espresso_settings_class = from_dict(data_class=EspressoSettings, 
                                        data=input_settings_dict,
                                        config=Config(type_hooks={str: str.lower}, strict=True))
    

    from ase.io.espresso import read_fortran_namelist    
    
    # parse namelist section and extract remaining lines
    with open(espresso_settings_class.pwi_path, 'r') as file:
        dftprogram_settings_dict, card_lines = read_fortran_namelist(file)

    dftprogram_settings_dict.update(vars(espresso_settings_class))

    #parse ATOMIC_SPECIES, K_POINTS and HUBBARD
    hubbard_corrections = False
    for i, line in enumerate(card_lines):
        if('ATOMIC_SPECIES' in line.upper()):
            atomic_species_index = i
        if('K_POINTS' in line.upper()):
            k_points_index = i
        if('HUBBARD' in line.upper()):
            hubbard_index = i
            hubbard_corrections = True

    def end_of_card(card, line):

        cards_list = ['ATOMIC_SPECIES',
                      'ATOMIC_POSITIONS',
                      'K_POINTS',
                      'ADDITIONAL_K_POINTS',
                      'CELL_PARAMETERS',
                      'CONSTRAINTS',
                      'OCCUPATIONS',
                      'ATOMIC_VELOCITIES',
                      'ATOMIC_FORCES',
                      'SOLVENTS',
                      'HUBBARD']

        for other_card in cards_list:
            if other_card == card.upper(): continue
            if other_card in line.upper(): return True
        
        return False

    #ATOMIC_SPECIES
    dftprogram_settings_dict['pseudopotentials'] = {}
    i = atomic_species_index+1
    while i < len(card_lines):
        line = card_lines[i]
        
        if end_of_card(card='ATOMIC_SPECIES', line=line): break       
        
        element, mass, pseudo = line.split()
        dftprogram_settings_dict['pseudopotentials'].update({element : pseudo})
        i+=1

    #K_POINTS
    if 'gamma' in card_lines[k_points_index].split()[1].strip().lower():
        dftprogram_settings_dict['kpts'] = None
        dftprogram_settings_dict['koffset'] = None
    else:
        line = card_lines[k_points_index+1]
        dftprogram_settings_dict['kpts'] = list(map(int, line.split()[:3]))
        dftprogram_settings_dict['koffset'] = list(map(int, line.split()[3:]))


    dftprogram_settings_dict['additional_cards'] = []

    #HUBBARD
    if hubbard_corrections:
        i = hubbard_index
        while i < len(card_lines):
            line = card_lines[i]
            
            if end_of_card(card='HUBBARD', line=line): break       
            
            dftprogram_settings_dict['additional_cards'] += [line]    
            i+=1

    return dftprogram_settings_dict




@dataclass
class VaspSettings:
    vasp_pp_path: str
    vasp_pseudo_setups: Optional[dict]
    pymatgen_set: Optional[str] 
    incar_path: Optional[str]
    kpoints_path: Optional[str]
    ediffg_screening: float = HYBRID_SCREENING_THRESHOLDS['VASP']   
    vasp_xc_functional: str = "PBE"   

def read_vasp_settings(input_settings_dict : dict):

    vasp_settings_class = from_dict(data_class=VaspSettings, 
                                    data=input_settings_dict,
                                    config=Config(type_hooks={str: str.lower}, strict=True))
    
    
    dftprogram_settings_dict = vars(vasp_settings_class)

    if vasp_settings_class.incar_path:
        with open(vasp_settings_class.incar_path, 'r') as f:
            dftprogram_settings_dict["incar_string"] = f.read()

    if vasp_settings_class.kpoints_path:
        with open(vasp_settings_class.kpoints_path, 'r') as f:
            dftprogram_settings_dict["kpoints_string"] = f.read()

    return dftprogram_settings_dict
    

def get_dftprogram_settings(program : str, input_settings_dict : dict):
    
    if program == 'ESPRESSO':
        return read_espresso_settings(input_settings_dict)
    elif program == 'VASP':
        return read_vasp_settings(input_settings_dict)
    else:
        raise ValueError(f'Program {program} not recognized.')
