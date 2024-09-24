
from xsorb.settings import Settings

def override_settings(settings : Settings, calc_type : str):
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
            settings.dftprogram_settings_dict['ions'].update({'upscale': 1}) #TODO: lasciarlo scegliere all'utente


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

def override_settings_isolated_fragment(settings : Settings, natoms_mol : int, manual_dft_override : dict = None):
    '''
    Overrides the DFT code settings with specific options for the isolated fragment calculation

    Args:
    - settings: Settings object, containing the dft code parameters
    - natoms_mol: number of atoms in the molecule
    - manual_dft_override: additional settings fragment-specific, read from fragments.json
    '''
    if settings.program == 'ESPRESSO':

        settings.dftprogram_settings_dict['control'].update({'calculation' : 'relax'})
        settings.dftprogram_settings_dict['control'].update({'restart_mode' : 'from_scratch'})
        settings.dftprogram_settings_dict['control'].update({'outdir' : 'WORK'}) #TODO: lasciarlo scegliere all'utente
        settings.dftprogram_settings_dict['system'].update({'nosym' : True})

        if settings.dftprogram_settings_dict['system'].get('nspin') == 2:
            for i, pseudo in enumerate(settings.dftprogram_settings_dict['pseudopotentials']):
                settings.dftprogram_settings_dict['system'].update({f'starting_magnetization({i+1})' : 0.6})        
        settings.dftprogram_settings_dict['electrons'].update({'mixing_beta' : 0.1})
        settings.dftprogram_settings_dict['kpts'] = None
        settings.dftprogram_settings_dict['koffset'] = None

        if manual_dft_override is not None:
    
            if 'SYSTEM' in manual_dft_override:
                if 'nspin' not in manual_dft_override['SYSTEM']:
                    settings.dftprogram_settings_dict['system'].update({'nspin' : 2}) #spin-polarized for isolated fragments (if not explicitly disabled in manual_dft_override)
                
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
            
            key = line.split('=')[0].strip()
            if key in manual_dft_override:
                s[i] = f'{key} = {manual_dft_override[key].strip()}'
                manual_dft_override.pop(key)
                
        if missing_isym: s.append('ISYM = 0')
        if missing_magmom: s.append(f'MAGMOM = {natoms_mol}*0.6')

        #TODO: check that if we set MAGMOM when ispin = 1 it does not give an error in VASP


        if manual_dft_override is not None:
            for key in manual_dft_override:
                s.append(f'{key} = {manual_dft_override[key].strip()}')


        settings.dftprogram_settings_dict['incar_string'] = '\n'.join(s)     

        settings.dftprogram_settings_dict['kpoints_string'] = 'Isolated fragment\n0\nGamma\n1 1 1\n'

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
            
            if card and dft_settings_dict and card in dft_settings_dict:
                found = False
                for key, val in dft_settings_dict[card].items():
                    if key in line: 
                        dft_newlines.append(f'   {key} = {val}\n')
                        dft_settings_dict[card].pop(key)
                        found = True
                        break
                if found: continue
                
            if card and dft_settings_dict and line.strip()=='/': 
                for key, val in dft_settings_dict[card].items():
                    dft_newlines.append(f'   {key} = {val}\n')
                dft_newlines.append(line)
                card = None

            dft_newlines.append(line)
            


    elif program == 'VASP':
        
        in_incar = False
        for line in dft_section:
            
            if '&INCAR' in line: 
                in_incar = True
                dft_newlines.append(line)
                continue
            
            if in_incar and dft_settings_dict :
                found = False
                for key, val in dft_settings_dict.items():
                    if key in line: 
                        dft_newlines.append(f'   {key} = {val}\n')
                        dft_settings_dict.pop(key)
                        found = True
                        break
                if found: continue
                
            if in_incar and dft_settings_dict  and line.strip()=='/': 
                for key, val in dft_settings_dict.items():
                    dft_newlines.append(f'   {key} = {val}\n')
                dft_newlines.append(line)
                card = None

            dft_newlines.append(line)

        #TODO: set the magmoms if not in dft_settings_dict, to match the number of atoms
        

    return dft_newlines    