import numpy as np
from ReferenceFiles.unit_conversions import *

#This is the data obtained from quantum esspresso
#pbe-nl-kjpaw_psl.1.0.0

# SCF cycle convergence threshold
convergence_threshold = 1     # in Ry = 1 meV

# Simulation details
chemical_formula = 'Si'
structure_names = ['$Fd\overline{3}m$', '$I4_{1}/amd$']  # Hermann-Mauguin notation for diamond & beta-Sn structures
n_atom = 16


volumes_sim_diamond = cubic_meters_per_cubic_angstrom*np.array([22.75991690336,22.99949497764,23.23907304949,23.47865112133,23.71822919318,23.95780726746,24.19738534176,24.43696341360,24.67654148545,24.91611955729,25.15569763157])
total_energies_strain_diamond = (joules_per_rydberg/n_atom)*np.array([-747.25298904,-747.26072309,-747.26668142,-747.27085720,-747.27337511,-747.27420258,-747.27342032,-747.27099990,-747.26703600,-747.26163010,-747.25471078])
volumes_sim_betasn = cubic_meters_per_cubic_angstrom*np.array([12.75576087798,14.57801243116,16.40026398432,18.22251553943,20.04476709452,21.86701864770,23.68927020279])
total_energies_strain_betasn = (joules_per_rydberg/n_atom)*np.array([-744.04950150,-746.23214087,-746.97811086,-747.13420595,-747.06624627,-746.94247460,-746.81741950])