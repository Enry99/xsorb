from pathlib import Path

p = Path('LICENSE')

#print(p.is_file())


from dataclasses import dataclass, field
from typing import Optional
@dataclass
class VaspSettings:
    vasp_pp_path: str
    ediffg_screening: float =  10  
    vasp_xc_functional: str = "PBE"   

vs = VaspSettings(vasp_pp_path = 'path/to/vasp_pp')

print(vars(vs))