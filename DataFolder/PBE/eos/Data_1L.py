import numpy as np
from ReferenceFiles.unit_conversions import *

#This is the data obtained from quantum esspresso
#pbe-nl-kjpaw_psl.1.0.0

# SCF cycle convergence threshold
convergence_threshold = 1     # in Ry = 1 meV

# Simulation details
chemical_formula = 'Si'
structure_names = ['$Fd\overline{3}m$', '$I4_{1}/amd$']  # Hermann-Mauguin notation for diamond & beta-Sn structures
n_atom = 8


# volumes_sim_diamond = cubic_meters_per_cubic_angstrom * np.array([20.12663366986, 20.33849296808, 20.55035227127, 20.76221156948, 20.97407087773, 21.18593017590, 21.39778947408, 21.60964878233, 21.82150808054, 22.03336738373, 22.24522668195])
# total_energies_strain_diamond = (joules_per_rydberg / n_atom) * np.array([-373.50358276, -373.50574327, -373.50735001, -373.50848663, -373.50932255, -373.50962104, -373.50929121, -373.50850353, -373.50725521, -373.50575021, -373.50360885])
volumes_sim_betasn = cubic_meters_per_cubic_angstrom * np.array([10.59331112283, 12.10664127782, 13.61997143660, 15.13330159916, 16.64663175794, 18.15996192050, 19.67329207550])
total_energies_strain_betasn = (joules_per_rydberg / n_atom) * np.array([-372.20766344, -373.15406152, -373.47631757, -373.54085602, -373.50807039, -373.44924960, -373.38985567])
volumes_sim_diamond = cubic_meters_per_cubic_angstrom*np.array([20.12926999308,20.34115704327,20.55304409848,20.76493114867,20.97681820388,21.18870525403,21.40059230421,21.61247935942,21.82436640961,22.03625346482,22.24814051497])
total_energies_strain_diamond = (joules_per_rydberg/n_atom)*np.array([-373.50366329,-373.50580660,-373.50748463,-373.50868629,-373.50940685,-373.50964594,-373.50940380,-373.50868015,-373.50747575,-373.50579797,-373.50367178])
