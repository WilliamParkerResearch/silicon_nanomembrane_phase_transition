import numpy as np
from ReferenceFiles.unit_conversions import *
n_atom = 48
cella_diamond = (1e-10)*np.array([5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790])
cellb_diamond = (1e-10)*np.array([5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790])
cellc_diamond = (1e-10)*np.array([5.21687452666666666666,5.27178899546666666666,5.32670346400000000000,5.38161793280000000000,5.43653240133333333333,5.49144687013333333333,5.54636133893333333333,5.60127580746666666666,5.65619027626666666666,5.71110474480000000000,5.76601921360000000000])
total_energies_strain_diamond = (joules_per_rydberg/n_atom)*np.array([-2242.32742408,-2242.35729830,-2242.38028488,-2242.39652435,-2242.40612181,-2242.40921959,-2242.40590161,-2242.39638899,-2242.38082314,-2242.35939825,-2242.33224812])
volumes_sim_diamond = cella_diamond*cellb_diamond*cellc_diamond/8
n_atom = 48
cella_betasn = (1e-10)*np.array([4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525])
cellb_betasn = (1e-10)*np.array([4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525])
cellc_betasn = (1e-10)*np.array([5.01966401006451612903,5.07250257858064516129,5.12534114709677419354,5.17817971561290322580,5.23101828412903225806,5.28385685264516129032,5.33669542116129032258,5.38953398967741935483,5.44237255819354838709,5.49521112670967741935,5.54804969522580645161])
total_energies_strain_betasn = (joules_per_rydberg/n_atom)*np.array([-2241.63129887,-2241.65885224,-2241.67886273,-2241.69235306,-2241.70008753,-2241.70253122,-2241.70012563,-2241.69316041,-2241.68210005,-2241.66748325,-2241.64985579])
volumes_sim_betasn = cella_betasn*cellb_betasn*cellc_betasn/8