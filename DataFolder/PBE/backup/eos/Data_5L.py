import numpy as np
from ReferenceFiles.unit_conversions import *

#This is the data obtained from quantum esspresso
#pbe-nl-kjpaw_psl.1.0.0

# SCF cycle convergence threshold
convergence_threshold = 1     # in Ry = 1 meV

# Simulation details
chemical_formula = 'Si'
structure_names = ['$Fd\overline{3}m$', '$I4_{1}/amd$']  # Hermann-Mauguin notation for diamond & beta-Sn structures
n_atom = 40

volumes_sim_diamond = cubic_meters_per_cubic_angstrom * np.array([24.27685529827, 24.53240114399, 24.78794698971, 25.04349283445, 25.29903867918, 25.55458452491, 25.81013037063, 26.06567621537, 26.32122206010, 26.57676790582, 26.83231375154])
total_energies_strain_diamond = (joules_per_rydberg / n_atom) * np.array([-1868.52741765, -1868.55180741, -1868.57055549, -1868.58381040, -1868.59165491, -1868.59412512, -1868.59146848, -1868.58370430, -1868.57101034, -1868.55355050, -1868.53140716])
volumes_sim_betasn = cubic_meters_per_cubic_angstrom * np.array([18.54036898650, 18.73553076484, 18.93069254395, 19.12585432306, 19.32101610141, 19.51617788052, 19.71133965963, 19.90650143798, 20.10166321709, 20.29682499620, 20.49198677455])
total_energies_strain_betasn = (joules_per_rydberg / n_atom) * np.array([-1867.99770172, -1868.02053312, -1868.03732584, -1868.04870033, -1868.05518266, -1868.05723175, -1868.05515575, -1868.04925380, -1868.03997991, -1868.02787716, -1868.01343270])