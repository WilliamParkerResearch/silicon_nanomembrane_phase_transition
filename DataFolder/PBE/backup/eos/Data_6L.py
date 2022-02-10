import numpy as np
from ReferenceFiles.unit_conversions import *

#This is the data obtained from quantum esspresso
#pbe-nl-kjpaw_psl.1.0.0

# SCF cycle convergence threshold
convergence_threshold = 1     # in Ry = 1 meV

# Simulation details
chemical_formula = 'Si'
structure_names = ['$Fd\overline{3}m$', '$I4_{1}/amd$']  # Hermann-Mauguin notation for diamond & beta-Sn structures
n_atom = 48


volumes_sim_diamond = cubic_meters_per_cubic_angstrom * np.array([24.42244082389, 24.67951914809, 24.93659747312, 25.19367579732, 25.45075412235, 25.70783244655, 25.96491077075, 26.22198909577, 26.47906741998, 26.73614574500, 26.99322406920])
total_energies_strain_diamond = (joules_per_rydberg / n_atom) * np.array([-2242.32742408, -2242.35729830, -2242.38028488, -2242.39652435, -2242.40612181, -2242.40921959, -2242.40590161, -2242.39638899, -2242.38082314, -2242.35939825, -2242.33224812])
volumes_sim_betasn = cubic_meters_per_cubic_angstrom * np.array([18.67629879901, 18.87289141824, 19.06948403747, 19.26607665543, 19.46266927466, 19.65926189389, 19.85585451313, 20.05244713236, 20.24903975032, 20.44563236955])
total_energies_strain_betasn = (joules_per_rydberg / n_atom) * np.array([-2241.63129887, -2241.65885224, -2241.67886273, -2241.69235306, -2241.70008753, -2241.70253122, -2241.70012563, -2241.69316041, -2241.68210005, -2241.66748325])