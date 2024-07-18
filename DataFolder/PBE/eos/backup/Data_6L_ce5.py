import numpy as np
from ReferenceFiles.unit_conversions import *
n_atom = 48
cella_diamond = (1e-10)*np.array([5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790])
cellb_diamond = (1e-10)*np.array([5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790])
cellc_diamond = (1e-10)*np.array([5.21797385867,5.27289989867,5.32782593867,5.38275198000,5.43767802000,5.49260406133,5.54753010267,5.60245614267,5.65738218400,5.71230822400,5.76723426400])
total_energies_strain_diamond = (joules_per_rydberg/n_atom)*np.array([-2242.32742408,-2242.35729830,-2242.38028488,-2242.39652435,-2242.40612181,-2242.40921959,-2242.40590161,-2242.39638899,-2242.38082314,-2242.35939825,-2242.33224812])
volumes_sim_diamond = cella_diamond*cellb_diamond*cellc_diamond/8
n_atom = 48
cella_betasn = (1e-10)*np.array([4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525])
cellb_betasn = (1e-10)*np.array([4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525])
cellc_betasn = (1e-10)*np.array([4.96380071771,5.01605125257,5.06830178629,5.12055231886,5.17280285371,5.22505338743,5.27730392114,5.32955445600,5.38180498857,5.43405552229,5.48630605714])
total_energies_strain_betasn = (joules_per_rydberg/n_atom)*np.array([-2241.63129887,-2241.65885224,-2241.67886273,-2241.69235306,-2241.70008753,-2241.70253122,-2241.70012563,-2241.69316041,-2241.68210005,-2241.66748325,-2241.64985579])
volumes_sim_betasn = cella_betasn*cellb_betasn*cellc_betasn/8