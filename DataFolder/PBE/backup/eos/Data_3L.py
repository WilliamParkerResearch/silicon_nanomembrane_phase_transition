import numpy as np
from ReferenceFiles.unit_conversions import *

#This is the data obtained from quantum esspresso
#pbe-nl-kjpaw_psl.1.0.0

# SCF cycle convergence threshold
convergence_threshold = 1     # in Ry = 1 meV

# Simulation details
chemical_formula = 'Si'
structure_names = ['$Fd\overline{3}m$', '$I4_{1}/amd$']  # Hermann-Mauguin notation for diamond & beta-Sn structures
n_atom = 24


volumes_sim_diamond = cubic_meters_per_cubic_angstrom * np.array([23.82678988521, 24.07759820030, 24.32840651540, 24.57921482886, 24.83002314397, 25.08083145906, 25.33163977416, 25.58244808927, 25.83325640272, 26.08406471782, 26.33487303291])
total_energies_strain_diamond = (joules_per_rydberg / n_atom) * np.array([-1121.05221904, -1121.06512793, -1121.07507526, -1121.08211431, -1121.08627268, -1121.08765665, -1121.08628433, -1121.08222083, -1121.07553173, -1121.06628301, -1121.05453921])
volumes_sim_betasn = cubic_meters_per_cubic_angstrom * np.array([18.13183567901, 18.32070896802, 18.50958225574, 18.69845554346, 18.88732883246, 19.07620212147, 19.26507540919, 19.45394869691, 19.64282198591])
total_energies_strain_betasn = (joules_per_rydberg / n_atom) * np.array([-1120.74629933, -1120.75651057, -1120.76334662, -1120.76713734, -1120.76833481, -1120.76723292, -1120.76415930, -1120.75932850, -1120.75283412])