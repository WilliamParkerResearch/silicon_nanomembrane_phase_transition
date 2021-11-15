import numpy as np
from ReferenceFiles.unit_conversions import *
n_atom = 8
cella_diamond = (1e-10)*np.array([5.49477513503,5.49477513503,5.49477513503,5.49477513503,5.49477513503,5.49477513503,5.49477513503,5.49477513503,5.49477513503,5.49477513503,5.49477513503])
cellb_diamond = (1e-10)*np.array([5.49477513503,5.49477513503,5.49477513503,5.49477513503,5.49477513503,5.49477513503,5.49477513503,5.49477513503,5.49477513503,5.49477513503,5.49477513503])
cellc_diamond = (1e-10)*np.array([5.33357201333333333333,5.38971487600000000000,5.44585774000000000000,5.50200060266666666666,5.55814346666666666666,5.61428632933333333333,5.67042919200000000000,5.72657205600000000000,5.78271491866666666666,5.83885778266666666666,5.89500064533333333333])
total_energies_strain_diamond = (joules_per_rydberg/n_atom)*np.array([-373.50366329,-373.50580660,-373.50748463,-373.50868629,-373.50940685,-373.50964594,-373.50940380,-373.50868015,-373.50747575,-373.50579797,-373.50367178])
volumes_sim_diamond = cella_diamond*cellb_diamond*cellc_diamond/8
n_atom = 8
cella_betasn = (1e-10)*np.array([4.77033373765,4.77033373765,4.77033373765,4.77033373765,4.77033373765,4.77033373765,4.77033373765])
cellb_betasn = (1e-10)*np.array([4.77033373765,4.77033373765,4.77033373765,4.77033373765,4.77033373765,4.77033373765,4.77033373765])
cellc_betasn = (1e-10)*np.array([3.19210668342857142857,3.64812192228571428571,4.10413716228571428571,4.56015240342857142857,5.01616764342857142857,5.47218288457142857142,5.92819812342857142857])
total_energies_strain_betasn = (joules_per_rydberg/n_atom)*np.array([-372.20766344,-373.15406152,-373.47631757,-373.54085602,-373.50807039,-373.44924960,-373.38985567])
volumes_sim_betasn = cella_betasn*cellb_betasn*cellc_betasn/8
