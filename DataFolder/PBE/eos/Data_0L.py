import numpy as np
from ReferenceFiles.unit_conversions import eV_per_rydberg, meV_per_rydberg, joules_per_rydberg, meters_per_angstrom

#This is the data obtained from quantum esspresso
#pbe-nl-kjpaw_psl.1.0.0

# SCF cycle convergence threshold
convergence_threshold = 1     # in Ry = 1 meV

# Simulation details
chemical_formula = 'Si'
structure_names = ['$Fd\overline{3}m$', '$I4_{1}/amd$']  # Hermann-Mauguin notation for diamond & beta-Sn structures
n_atom = 8

# Diamond
total_energies_strain_diamond = np.array([-373.70545147, -373.74029300, -373.76659418, -373.78458837, -373.79529060, -373.79925178, -373.79734561, -373.78974896, -373.77718787, -373.76012605, -373.73918952])
lattice_parameters_diamond = np.array([9.81537, 9.91869, 10.022, 10.1253, 10.2286, 10.332, 10.4353, 10.5386, 10.6419, 10.7452, 10.8486])
celldm_1_diamond = 10.33196444
celldm_3_diamond = 1

#(betasn
total_energies_strain_all_betasn = np.array([-186.78113214, -186.79695038, -186.80862965, -186.81655692, -186.82108286, -186.82253892, -186.82122818, -186.81744443, -186.81146125, -186.80352771, -186.79387086])
lattice_parameters_all_betasn = np.array([8.60693, 8.69753, 8.78813, 8.87873, 8.96933, 9.05993, 9.15053, 9.24113, 9.33173, 9.42232, 9.51292])
total_energies_strain_shape_betasn = np.array([-186.78174938, -186.79733717, -186.80884250, -186.81664050, -186.82110034, -186.82253982, -186.82126390, -186.81755527, -186.81168132, -186.80389893, -186.79449005])
lattice_parameters_shape_betasn = np.array([8.60693, 8.69753, 8.78813, 8.87873, 8.96933, 9.05993, 9.15053, 9.24113, 9.33173, 9.42232, 9.51292])
celldm_1_betasn = 9.059927516
celldm_3_betasn = 0.556964722

# Convert array units
total_energies_strain_diamond = joules_per_rydberg * total_energies_strain_diamond / n_atom         # J / atom
total_energies_strain_betasn = joules_per_rydberg * total_energies_strain_all_betasn / (n_atom/2)


lattice_parameters_diamond = meters_per_angstrom * lattice_parameters_diamond                               # m
lattice_parameters_betasn = meters_per_angstrom * lattice_parameters_all_betasn
celldm_1_diamond = meters_per_angstrom * celldm_1_diamond
celldm_1_betasn = meters_per_angstrom * celldm_1_betasn
volumes_sim_diamond = celldm_3_diamond*(np.power(lattice_parameters_diamond,3))/n_atom
volumes_sim_betasn = celldm_3_betasn*(np.power(lattice_parameters_betasn,3))/(n_atom/2)
print(volumes_sim_betasn)
print(volumes_sim_diamond)