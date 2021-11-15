import numpy as np
from ReferenceFiles.unit_conversions import *
n_atom = 32
cella_diamond = (1e-10)*np.array([5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790])
cellb_diamond = (1e-10)*np.array([5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790,5.44010075790])
cellc_diamond = (1e-10)*np.array([5.22158766133,5.27655174267,5.33151582400,5.38647990533,5.44144398400,5.49640806533,5.55137214667,5.60633622533,5.66130030667,5.71626438800,5.77122846933])
total_energies_strain_diamond = (joules_per_rydberg/n_atom)*np.array([-1494.72712190,-1494.74605927,-1494.76059731,-1494.77090485,-1494.77702169,-1494.77900039,-1494.77691928,-1494.77094116,-1494.76113125,-1494.74767171,-1494.73060656])
volumes_sim_diamond = cella_diamond*cellb_diamond*cellc_diamond/8
n_atom = 32
cella_betasn = (1e-10)*np.array([4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525])
cellb_betasn = (1e-10)*np.array([4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525])
cellc_betasn = (1e-10)*np.array([4.95004270057,5.00214841257,5.05425412571,5.10635983886,5.15846555086,5.21057126400,5.26267697714,5.31478268914,5.36688840229,5.41899411543,5.47109982743])
total_energies_strain_betasn = (joules_per_rydberg/n_atom)*np.array([-1494.36684143,-1494.38435797,-1494.39698514,-1494.40545457,-1494.41032818,-1494.41189467,-1494.41037109,-1494.40588925,-1494.39870156,-1494.38914426,-1494.37759102])
volumes_sim_betasn = cella_betasn*cellb_betasn*cellc_betasn/8
