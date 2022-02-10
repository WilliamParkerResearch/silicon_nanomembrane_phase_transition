import numpy as np
from ReferenceFiles.unit_conversions import *

#This is the data obtained from quantum esspresso
#pbe-nl-kjpaw_psl.1.0.0

# SCF cycle convergence threshold
convergence_threshold = 1     # in Ry = 1 meV

# Simulation details
chemical_formula = 'Si'
structure_names = ['$Fd\overline{3}m$', '$I4_{1}/amd$']  # Hermann-Mauguin notation for diamond & beta-Sn structures
n_atom = 32


volumes_sim_diamond = cubic_meters_per_cubic_angstrom * np.array([23.91483998890, 24.16657514730, 24.41831030570, 24.67004546410, 24.92178062003, 25.17351577844, 25.42525093684, 25.67698609277, 25.92872125117, 26.18045640957, 26.43219156798])
total_energies_strain_diamond = (joules_per_rydberg / n_atom) * np.array([-1494.72712190, -1494.74605927, -1494.76059731, -1494.77090485, -1494.77702169, -1494.77900039, -1494.77691928, -1494.77094116, -1494.76113125, -1494.74767171, -1494.73060656])
volumes_sim_betasn = cubic_meters_per_cubic_angstrom*np.array([18.35159522610,18.54476991291,18.73794459974,18.93111928560,19.12429397243,19.31746865925,19.51064334607,19.70381803289,19.89699271875,20.09016740558])
total_energies_strain_betasn = (joules_per_rydberg/n_atom)*np.array([-1494.36684143,-1494.38435797,-1494.39698514,-1494.40545457,-1494.41032818,-1494.41189467,-1494.41037109,-1494.40588925,-1494.39870156,-1494.38914426])