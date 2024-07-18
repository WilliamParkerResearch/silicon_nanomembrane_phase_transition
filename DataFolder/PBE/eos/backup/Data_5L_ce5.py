import numpy as np
from ReferenceFiles.unit_conversions import *
n_atom = 40
cella_diamond = (1e-10)*np.array([5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790])
cellb_diamond = (1e-10)*np.array([5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790])
cellc_diamond = (1e-10)*np.array([5.36504622400,5.42152039600,5.47799456667,5.53446873733,5.59094290800,5.64741707867,5.70389124933,5.76036542000,5.81683959067,5.87331376133,5.92978793333])
total_energies_strain_diamond = (joules_per_rydberg/n_atom)*np.array([-1868.52741765,-1868.55180741,-1868.57055549,-1868.58381040,-1868.59165491,-1868.59412512,-1868.59146848,-1868.58370430,-1868.57101034,-1868.55355050,-1868.53140716])
volumes_sim_diamond = cella_diamond*cellb_diamond*cellc_diamond/8
n_atom = 40
cella_betasn = (1e-10)*np.array([4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525])
cellb_betasn = (1e-10)*np.array([4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525])
cellc_betasn = (1e-10)*np.array([4.74952804229,4.79952307429,4.84951810629,4.89951313829,4.94950817029,4.99950320229,5.04949823429,5.09949326629,5.14948829829,5.19948333029,5.24947836229])
total_energies_strain_betasn = (joules_per_rydberg/n_atom)*np.array([-1867.99770172,-1868.02053312,-1868.03732584,-1868.04870033,-1868.05518266,-1868.05723175,-1868.05515575,-1868.04925380,-1868.03997991,-1868.02787716,-1868.01343270])
volumes_sim_betasn = cella_betasn*cellb_betasn*cellc_betasn/8