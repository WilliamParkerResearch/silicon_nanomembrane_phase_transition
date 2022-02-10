import numpy as np
import matplotlib.pyplot as plt
from eos_information import eos_properties
from ReferenceFiles.cohesive_energy_values import *
from ReferenceFiles.unit_conversions import *
from scipy.optimize import curve_fit, fmin
from ReferenceFiles.FunctionDefinitions import mid, square_differences
from ReferenceFiles.plot_formating import *
from matplotlib.lines import Line2D

XC = 'SCAN'
fold = 'eosz'
# showplot=True
showplot=False
prd = 1
prb = 1
# #Use TeX fonts
# mpl.rcParams['text.usetex'] = True
# mpl.rcParams['font.sans-serif'] = "cmr10"
# mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['axes.formatter.useoffset'] = False

adj_l = 1.2
linewidth = universal_linewidth

n_layers = np.arange(0,5,1)
# n_layers = np.array([0,1,3,4])
# n_layers = np.array([0,1,2,7])

transition_pressures = np.array([])
transition_volumes_diamond = np.array([])
transition_volumes_betasn = np.array([])
volumes0_diamond = np.array([])
volumes0_betasn = np.array([])
energy0_diamond = np.array([])
energy0_betasn = np.array([])
bulk_mudulus0_diamond = np.array([])
bulk_mudulus0_betasn = np.array([])
cell_a_diamond = np.array([])
cell_c_diamond = np.array([])
cell_a_betasn = np.array([])
cell_c_betasn = np.array([])

for i in n_layers:
    print(i," layer")
    if i == 0:
        t_pres0, vol_d0, vol_b0, fitp_d0, fitp_b0,ad0,ab0,cd0,cb0,_,_ = eos_properties(i,exchange_correlation=XC,eos_folder=fold,propd=prd,propb=prb)
        vol0_d0 = fitp_d0[3]
        energy0_d0 = fitp_d0[0]
        bulk_mod0_d0 = fitp_d0[1]
        vol0_b0 = fitp_b0[3]
        energy0_b0 = fitp_b0[0]
        bulk_mod0_b0 = fitp_b0[1]
    else:
        t_pres, vol_d, vol_b, fitp_d, fitp_b,ad,ab,cd,cb,_,_ = eos_properties(i,exchange_correlation=XC,eos_folder=fold,propd=prd,propb=prb)
        transition_pressures = np.append(transition_pressures,t_pres)
        transition_volumes_diamond = np.append(transition_volumes_diamond,vol_d)
        transition_volumes_betasn = np.append(transition_volumes_betasn,vol_b)
        volumes0_diamond = np.append(volumes0_diamond,fitp_d[3])
        volumes0_betasn = np.append(volumes0_betasn,fitp_b[3])
        energy0_diamond = np.append(energy0_diamond,fitp_d[0])
        energy0_betasn = np.append(energy0_betasn,fitp_b[0])
        bulk_mudulus0_diamond = np.append(bulk_mudulus0_diamond,fitp_d[1])
        bulk_mudulus0_betasn = np.append(bulk_mudulus0_betasn,fitp_b[1])
        cell_a_diamond = np.append(cell_a_diamond,ad)
        cell_a_betasn = np.append(cell_a_betasn,ab)
        cell_c_diamond = np.append(cell_c_diamond, cd)
        cell_c_betasn = np.append(cell_c_betasn, cb)


if XC == 'PBE':
    atom_energy = atom_energy_pbe
elif XC == 'SCAN':
    atom_energy = atom_energy_scan

cohesive_energy0_diamond = -((energy0_diamond)-atom_energy)
cohesive_energy0_d0 = -((energy0_d0)-atom_energy)
cohesive_energy0_betasn = -((energy0_betasn)-atom_energy)
cohesive_energy0_b0 = -((energy0_b0)-atom_energy)



pseudo_array = np.zeros((5,13))

N_MLT = 5
n_layerst = np.arange(N_MLT+1)
for i in n_layerst:
    if i == 0:
        pseudo_array[0][1] = round(vol0_d0*1e30,2)
        pseudo_array[1][1] = round(vol_d0*1e30,2)
        pseudo_array[2][1] = round(bulk_mod0_d0*1e-9,2)
        pseudo_array[3][1] = round(cohesive_energy0_d0*eV_per_joule,2)
        pseudo_array[4][1] = round(-t_pres0*1e-9,2)
    elif i > (n_layers[-1]):
        pseudo_array[0][i+1] = 0
        pseudo_array[1][i+1] = 0
        pseudo_array[2][i+1] = 0
        pseudo_array[3][i+1] = 0
        pseudo_array[4][i+1] = 0
    else:
        pseudo_array[0][i+1] = round(volumes0_diamond[i-1]*1e30,2)
        pseudo_array[1][i+1] = round(transition_volumes_diamond[i-1]*1e30,2)
        pseudo_array[2][i+1] = round(bulk_mudulus0_diamond[i-1]*1e-9,2)
        pseudo_array[3][i+1] = round(cohesive_energy0_diamond[i-1]*eV_per_joule,2)
        pseudo_array[4][i+1] = round(-transition_pressures[i-1]*1e-9,2)



for i in n_layerst:
    if i == 0:
        pseudo_array[0][7] = round(vol0_b0*1e30,2)
        pseudo_array[1][7] = round(vol_b0*1e30,2)
        pseudo_array[2][7] = round(bulk_mod0_b0*1e-9,2)
        pseudo_array[3][7] = round(cohesive_energy0_b0*eV_per_joule,2)
        pseudo_array[4][7] = round(-t_pres0*1e-9,2)
    elif i > (n_layers[-1]):
        pseudo_array[0][i+7] = 0
        pseudo_array[1][i+7] = 0
        pseudo_array[2][i+7] = 0
        pseudo_array[3][i+7] = 0
        pseudo_array[4][i+7] = 0
    else:
        pseudo_array[0][i+7] = round(volumes0_betasn[i-1]*1e30,2)
        pseudo_array[1][i+7] = round(transition_volumes_betasn[i-1]*1e30,2)
        pseudo_array[2][i+7] = round(bulk_mudulus0_betasn[i-1]*1e-9,2)
        pseudo_array[3][i+7] = round(cohesive_energy0_betasn[i-1]*eV_per_joule,2)
        pseudo_array[4][i+7] = round(-transition_pressures[i-1]*1e-9,2)

pseudo_array_str = pseudo_array.astype(str)


pseudo_array_str[0][0] = '& $V_0$ ($\AA^3$/atom)'
pseudo_array_str[1][0] = '& $V_{\\rm T}$ ($\AA^3$/atom)'
pseudo_array_str[2][0] = '& $K_{0}$ (GPa)'
pseudo_array_str[3][0] = '& $E_{\\rm coh}$ (eV/atom)'
pseudo_array_str[4][0] = '& $P_{\\rm T}$ (GPa)'

empty_idx = np.argwhere(pseudo_array_str == '0.0')
for i in empty_idx:
    print(i[0])
    pseudo_array_str[i[0]][i[1]] = ' '


row_idx = np.arange(5)

for i in row_idx:
    if i < row_idx[-1]:
        print("     & ".join(pseudo_array_str[i])+'    \\\\ \cline{2-14}')
    elif i == row_idx[-1]:
        print("     & ".join(pseudo_array_str[i])+'     \\\\ \hline\hline')
