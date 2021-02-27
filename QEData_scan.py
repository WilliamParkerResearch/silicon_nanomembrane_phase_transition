import numpy as np
#from unit_conversions import eV_per_rydberg, meV_per_rydberg, joules_per_rydberg, meters_per_angstrom
from unit_conversions import *

#This is the data obtained from quantum esspresso
#SCAN

# SCF cycle convergence threshold
convergence_threshold = 1     # in Ry = 1 meV

# Simulation details
chemical_formula = 'Si'
structure_names = ['$Fd\overline{3}m$', '$I4_{1}/amd$']  # Hermann-Mauguin notation for diamond & beta-Sn structures
n_atom_diamond = 8
n_atom_betasn = 4



# Diamond
total_energies_ecut_diamond = np.array([-62.73361926, -62.96709272, -62.98735791, -62.99860656, -63.00004283, -63.00052766, -63.00130887, -63.00171886, -63.00187468])
cutoff_energies_diamond = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90])
total_energies_kpoint_diamond = np.array([-62.91266767, -62.96114937, -62.96709272, -62.96798439, -62.96820946, -62.96817902, -62.96822198, -62.96823581, -62.96821538, -62.96820512, -62.96819526, -62.96819736, -62.96819049, -62.96818873, -62.96820338])
kpoints_diamond = np.array([2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16])
total_energies_strain_diamond = np.array([-62.86693326, -62.90444191, -62.93231053, -62.95169590, -62.96291019, -62.96709272, -62.90149084])
lattice_parameters_diamond = np.array([9.70335, 9.80549, 9.90763, 10.0098, 10.1119, 10.2141, 10.7248])
celldm_1_diamond = 10.21405661
celldm_3_diamond = 1
                ##density of state
dos_energies_diamond = np.array([])
density_diamond = np.array([])
fermi_energy_diamond = 6.3010
                ##bands
n_bands_diamond = 1
ticks_diamond = []
ticklabels_diamond = ["$\Gamma$","X","M","$\Gamma$","R","X","M","R"]


kpoint_indexes_diamond = np.array([])
bands_eigenvalues_diamond_1 = np.array([])




#BetaSn
total_energies_ecut_BetaSn = np.array([-31.21574121, -31.34306546, -31.35357862, -31.35249278, -31.36128360, -31.36179966, -31.36199229, -31.36213487, -31.36204074])
cutoff_energies_BetaSn = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90])
total_energies_kpoint_BetaSn = np.array([-31.47184020, -31.38453166, -31.35698483, -31.34802414, -31.34990701, -31.36467148, -31.36129581, -31.36258376, -31.36081379, -31.36003274, -31.36094660, -31.36075722, -31.36130035, -31.36128360, -31.36093534])
kpoints_BetaSn = np.array([2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16])
total_energies_strain_all_BetaSn = np.array([-31.31158251, -31.32813412, -31.34034594, -31.34858344, -31.35324402, -31.36128360, -31.35993836, -31.35607757, -31.34993719, -31.34176200, -31.33197332])
lattice_parameters_all_BetaSn = np.array([8.47207, 8.56124, 8.65042, 8.7396, 8.82878, 8.91796, 9.00714, 9.09632, 9.1855, 9.27468, 9.36386])
total_energies_strain_shape_BetaSn = np.array([-31.31293974, -31.32893368, -31.34075288, -31.34874031, -31.35327325, -31.36128896, -31.36003545, -31.35641082, -31.35079379, -31.34317682, -31.33481494])
lattice_parameters_shape_BetaSn = np.array([8.47207, 8.56124, 8.65042, 8.7396, 8.82878, 8.91796, 9.00714, 9.09632, 9.1855, 9.27468, 9.36386])
celldm_1_BetaSn = 8.917963303
celldm_3_BetaSn = 0.571032641
                ##density of state
dos_energies_betasn = np.array([])
density_betasn = np.array([])
fermi_energy_betasn = 9.8707
                ##bands
n_bands_betasn = 1
ticks_betasn = []
ticklabels_betasn = ["$\Gamma$","X","M","$\Gamma$","Z","R","A","Z","X","R","M","A"]


kpoint_indexes_betasn = np.array([])
bands_eigenvalues_betasn_1 = np.array([])

# Convert array units
total_energies_ecut_diamond = meV_per_rydberg * total_energies_ecut_diamond / n_atom_diamond  # meV / atom
total_energies_kpoint_diamond = meV_per_rydberg * total_energies_kpoint_diamond / n_atom_diamond
total_energies_ecut_BetaSn = meV_per_rydberg * total_energies_ecut_BetaSn / n_atom_betasn
total_energies_kpoint_BetaSn = meV_per_rydberg * total_energies_kpoint_BetaSn / n_atom_betasn
total_energies_strain_diamond = joules_per_rydberg * total_energies_strain_diamond / n_atom_diamond  # J / atom
total_energies_strain_all_BetaSn = joules_per_rydberg * total_energies_strain_all_BetaSn / n_atom_betasn
total_energies_strain_shape_BetaSn = joules_per_rydberg * total_energies_strain_shape_BetaSn / n_atom_betasn
cutoff_energies_diamond = eV_per_rydberg * cutoff_energies_diamond  # eV
cutoff_energies_BetaSn = eV_per_rydberg * cutoff_energies_BetaSn
lattice_parameters_diamond = meters_per_angstrom * lattice_parameters_diamond  # m
lattice_parameters_BetaSn = meters_per_angstrom * lattice_parameters_all_BetaSn
celldm_1_diamond = meters_per_angstrom * celldm_1_diamond
celldm_1_BetaSn = meters_per_angstrom * celldm_1_BetaSn
