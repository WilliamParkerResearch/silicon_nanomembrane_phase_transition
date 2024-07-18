import numpy as np
from ReferenceFiles.unit_conversions import *
n_atom = 64
cella_diamond = (1e-10)*np.array([5.46197051087,5.46197051087,5.46197051087,5.46197051087,5.46197051087,5.46197051087,5.46197051087,5.46197051087,5.46197051087,5.46197051087,5.46197051087])
cellb_diamond = (1e-10)*np.array([5.46197051087,5.46197051087,5.46197051087,5.46197051087,5.46197051087,5.46197051087,5.46197051087,5.46197051087,5.46197051087,5.46197051087,5.46197051087])
cellc_diamond = (1e-10)*np.array([36.146397314,38.048839278,39.634207581,40.860007815,41.693885525,42.115035884,42.536186243,43.386909968,44.688517267,46.476057957,48.799860855])
total_energies_strain_diamond = (joules_per_rydberg/n_atom)*np.array([-2989.14108060,-2989.62083046,-2989.88569163,-2990.00138510,-2990.03615707,-2990.04061689,-2990.03650593,-2990.00285034,-2989.88916121,-2989.62233769,-2989.14111239])
volumes_sim_diamond = cella_diamond*cellb_diamond*cellc_diamond/n_atom
n_atom = 64
cella_betasn = (1e-10)*np.array([4.80869776838,4.80869776838,4.80869776838,4.80869776838,4.80869776838,4.80869776838,4.80869776838,4.80869776838,4.80869776838,4.80869776838])
cellb_betasn = (1e-10)*np.array([4.80869776838,4.80869776838,4.80869776838,4.80869776838,4.80869776838,4.80869776838,4.80869776838,4.80869776838,4.80869776838,4.80869776838])
cellc_betasn = (1e-10)*np.array([37.312573492,38.867264055,40.069344387,40.887086110,41.300086979,41.713087848,42.547349605,43.823770093,45.576720897,47.855556942])
total_energies_strain_betasn = (joules_per_rydberg/n_atom)*np.array([-2988.53419726,-2988.85124400,-2988.96256763,-2988.99150100,-2988.99462592,-2988.99101245,-2988.96571966,-2988.89085658,-2988.74568708,-2988.51461409])
volumes_sim_betasn = cella_betasn*cellb_betasn*cellc_betasn/n_atom