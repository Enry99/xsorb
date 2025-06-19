'''
General settings and rendering parameters class for the atomsplot package.
'''

import os
import json
import logging
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from ase.data.colors import jmol_colors, cpk_colors

logger = logging.getLogger(__name__)


ATOMIC_RADIUS_DEFAULT = 0.7
BOND_RADIUS_DEFAULT = 0.8
BOND_LINE_WIDTH_DEFAULT = 0.15
CELL_LINE_WIDTH_DEFAULT = 0.05

vesta_colors = np.array([
    [1.00000, 0.00000, 0.00000],
    [1.00000, 0.80000, 0.80000],
    [0.98907, 0.91312, 0.81091],
    [0.52731, 0.87953, 0.45670],
    [0.37147, 0.84590, 0.48292],
    [0.12490, 0.63612, 0.05948],
    [0.50430, 0.28659, 0.16236],
    [0.69139, 0.72934, 0.90280],
    [0.99997, 0.01328, 0.00000],
    [0.69139, 0.72934, 0.90280],
    [0.99954, 0.21788, 0.71035],
    [0.97955, 0.86618, 0.23787],
    [0.98773, 0.48452, 0.08470],
    [0.50718, 0.70056, 0.84062],
    [0.10596, 0.23226, 0.98096],
    [0.75557, 0.61256, 0.76425],
    [1.00000, 0.98071, 0.00000],
    [0.19583, 0.98828, 0.01167],
    [0.81349, 0.99731, 0.77075],
    [0.63255, 0.13281, 0.96858],
    [0.35642, 0.58863, 0.74498],
    [0.71209, 0.38930, 0.67279],
    [0.47237, 0.79393, 1.00000],
    [0.90000, 0.10000, 0.00000],
    [0.00000, 0.00000, 0.62000],
    [0.66148, 0.03412, 0.62036],
    [0.71051, 0.44662, 0.00136],
    [0.00000, 0.00000, 0.68666],
    [0.72032, 0.73631, 0.74339],
    [0.13390, 0.28022, 0.86606],
    [0.56123, 0.56445, 0.50799],
    [0.62292, 0.89293, 0.45486],
    [0.49557, 0.43499, 0.65193],
    [0.45814, 0.81694, 0.34249],
    [0.60420, 0.93874, 0.06122],
    [0.49645, 0.19333, 0.01076],
    [0.98102, 0.75805, 0.95413],
    [1.00000, 0.00000, 0.60000],
    [0.00000, 1.00000, 0.15259],
    [0.40259, 0.59739, 0.55813],
    [0.00000, 1.00000, 0.00000],
    [0.29992, 0.70007, 0.46459],
    [0.70584, 0.52602, 0.68925],
    [0.80574, 0.68699, 0.79478],
    [0.81184, 0.72113, 0.68089],
    [0.80748, 0.82205, 0.67068],
    [0.75978, 0.76818, 0.72454],
    [0.72032, 0.73631, 0.74339],
    [0.95145, 0.12102, 0.86354],
    [0.84378, 0.50401, 0.73483],
    [0.60764, 0.56052, 0.72926],
    [0.84627, 0.51498, 0.31315],
    [0.67958, 0.63586, 0.32038],
    [0.55914, 0.12200, 0.54453],
    [0.60662, 0.63218, 0.97305],
    [0.05872, 0.99922, 0.72578],
    [0.11835, 0.93959, 0.17565],
    [0.35340, 0.77057, 0.28737],
    [0.82055, 0.99071, 0.02374],
    [0.99130, 0.88559, 0.02315],
    [0.98701, 0.55560, 0.02744],
    [0.00000, 0.00000, 0.96000],
    [0.99042, 0.02403, 0.49195],
    [0.98367, 0.03078, 0.83615],
    [0.75325, 0.01445, 1.00000],
    [0.44315, 0.01663, 0.99782],
    [0.19390, 0.02374, 0.99071],
    [0.02837, 0.25876, 0.98608],
    [0.28688, 0.45071, 0.23043],
    [0.00000, 0.00000, 0.88000],
    [0.15323, 0.99165, 0.95836],
    [0.15097, 0.99391, 0.71032],
    [0.70704, 0.70552, 0.35090],
    [0.71952, 0.60694, 0.33841],
    [0.55616, 0.54257, 0.50178],
    [0.70294, 0.69401, 0.55789],
    [0.78703, 0.69512, 0.47379],
    [0.78975, 0.81033, 0.45049],
    [0.79997, 0.77511, 0.75068],
    [0.99628, 0.70149, 0.22106],
    [0.82940, 0.72125, 0.79823],
    [0.58798, 0.53854, 0.42649],
    [0.32386, 0.32592, 0.35729],
    [0.82428, 0.18732, 0.97211],
    [0.00000, 0.00000, 1.00000],
    [0.00000, 0.00000, 1.00000],
    [1.00000, 1.00000, 0.00000],
    [0.00000, 0.00000, 0.00000],
    [0.42959, 0.66659, 0.34786],
    [0.39344, 0.62101, 0.45034],
    [0.14893, 0.99596, 0.47106],
    [0.16101, 0.98387, 0.20855],
    [0.47774, 0.63362, 0.66714],
    [0.30000, 0.30000, 0.30000],
    [0.30000, 0.30000, 0.30000],
    [0.30000, 0.30000, 0.30000]
])


@dataclass
class CustomSettings:
    """Configuration settings for rendering atoms
    """

    atomic_colors : dict = field(default_factory=dict)
    molecule_colors : dict = field(default_factory=dict)
    color_scheme : np.ndarray = field(default_factory=lambda: jmol_colors.copy())  # Use ASE's Jmol colors by default

    mol_indices : Optional[list[int]] = None  # Indices of molecules to highlight
    nontransparent_atoms : list[int] = field(default_factory=list)

    atomic_radius : float = ATOMIC_RADIUS_DEFAULT
    bond_radius : float = BOND_RADIUS_DEFAULT
    bond_line_width : float = BOND_LINE_WIDTH_DEFAULT
    cell_line_width : float = CELL_LINE_WIDTH_DEFAULT


    def __post_init__(self):
        """Post-initialization to set custom settings."""

        custom_settings_path = "image_settings.json"
        if os.path.isfile("image_settings.json"):
            with open(custom_settings_path, "r") as f:
                custom_settings : dict = json.load(f)
                logger.info("Custom colors read from %s file.", custom_settings_path)

            if custom_settings.pop("povray_old_style", None):
                os.environ['POVRAY_OLD_STYLE'] = '1'

            color_scheme = custom_settings.pop("color_scheme", "jmol").lower()
            if color_scheme == "vesta":
                self.color_scheme = vesta_colors.copy()
                logger.info("Using VESTA colors for elements.")
            elif color_scheme == "cpk":
                self.color_scheme = cpk_colors.copy()
                logger.info("Using CPK colors for elements.")
            elif color_scheme == "jmol":
                self.color_scheme = jmol_colors.copy()
                logger.info("Using Jmol colors for elements.")
            else:
                logger.warning("Unknown color scheme '%s'. Using Jmol colors by default.", color_scheme)

            for key, value in custom_settings.items():
                if hasattr(self, key):
                    setattr(self, key, value)
                else:
                    logger.warning("Custom setting '%s' not recognized.", key)
