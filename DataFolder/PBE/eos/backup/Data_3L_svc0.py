import numpy as np
from ReferenceFiles.unit_conversions import *
n_atom = 24
cella_diamond = (1e-10)*np.array([5.43923438402,5.43923438402,5.43923438402,5.43923438402,5.43923438402,5.43923438402,5.43923438402,5.43923438402,5.43923438402,5.43923438402,5.43923438402])
cellb_diamond = (1e-10)*np.array([5.43923438402,5.43923438402,5.43923438402,5.43923438402,5.43923438402,5.43923438402,5.43923438402,5.43923438402,5.43923438402,5.43923438402,5.43923438402])
cellc_diamond = (1e-10)*np.array([5.27144645636363636363,5.32693536654545454545,5.38242427672727272727,5.43791318654545454545,5.49340209672727272727,5.54889100690909090909,5.60437991709090909090,5.65986882727272727272,5.71535773709090909090,5.77084664727272727272,5.82633555745454545454])
total_energies_strain_diamond = (joules_per_rydberg/n_atom)*np.array([-1121.05221904,-1121.06512793,-1121.07507526,-1121.08211431,-1121.08627268,-1121.08765665,-1121.08628433,-1121.08222083,-1121.07553173,-1121.06628301,-1121.05453921])
volumes_sim_diamond = cella_diamond*cellb_diamond*cellc_diamond/8
n_atom = 24
cella_betasn = (1e-10)*np.array([4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525])
cellb_betasn = (1e-10)*np.array([4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525,4.79430526525])
cellc_betasn = (1e-10)*np.array([4.88739496347826086956,4.93884122608695652173,4.99028748904347826086,5.04173375165217391304,5.09318001426086956521,5.14462627721739130434,5.19607254017391304347,5.24751880278260869565,5.29896506539130434782,5.35041132834782608695,5.40185759095652173913])
total_energies_strain_betasn = (joules_per_rydberg/n_atom)*np.array([-1120.73225958,-1120.74629933,-1120.75651057,-1120.76334662,-1120.76713734,-1120.76833481,-1120.76723292,-1120.76415930,-1120.75932850,-1120.75283412,-1120.74488480])
volumes_sim_betasn = cella_betasn*cellb_betasn*cellc_betasn/8