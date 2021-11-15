import numpy as np
from ReferenceFiles.unit_conversions import *
n_atom = 40
cella_diamond = (1e-10)*np.array([5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790])
cellb_diamond = (1e-10)*np.array([5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790])
cellc_diamond = (1e-10)*np.array([5.18091167031578947368,5.23544758273684210526,5.28998349515789473684,5.34451940736842105263,5.39905531957894736842,5.45359123200000000000,5.50812714442105263157,5.56266305663157894736,5.61719896884210526315,5.67173488126315789473,5.72627079368421052631])
total_energies_strain_diamond = (joules_per_rydberg/n_atom)*np.array([-1868.52741765,-1868.55180741,-1868.57055549,-1868.58381040,-1868.59165491,-1868.59412512,-1868.59146848,-1868.58370430,-1868.57101034,-1868.55355050,-1868.53140716])
volumes_sim_diamond = cella_diamond*cellb_diamond*cellc_diamond/8
n_atom = 40
cella_betasn = (1e-10)*np.array([4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525])
cellb_betasn = (1e-10)*np.array([4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525])
cellc_betasn = (1e-10)*np.array([4.96379276656410256410,5.01604321661538461538,5.06829366687179487179,5.12054411712820512820,5.17279456717948717948,5.22504501743589743589,5.27729546769230769230,5.32954591774358974358,5.38179636800000000000,5.43404681825641025641,5.48629726830769230769])
total_energies_strain_betasn = (joules_per_rydberg/n_atom)*np.array([-1867.99770172,-1868.02053312,-1868.03732584,-1868.04870033,-1868.05518266,-1868.05723175,-1868.05515575,-1868.04925380,-1868.03997991,-1868.02787716,-1868.01343270])
volumes_sim_betasn = cella_betasn*cellb_betasn*cellc_betasn/8
