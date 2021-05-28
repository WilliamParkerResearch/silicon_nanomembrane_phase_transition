import scipy.optimize as sp
import numpy as np
import matplotlib.pyplot as plt
from ReferenceFiles.FunctionDefinitions import *
from ReferenceFiles.unit_conversions import *

# Parameters
number_of_volume_points = 100
number_of_enthalpy_samples = 10000
pressure_extreme = 50*11.7e10
exchange_correlation = 'PBE'
N_ML = '0'

directoryofdata='DataFolder'+'.'+exchange_correlation+'.'+'Data_'+N_ML+'L'
exec(f'from {directoryofdata} import *')

# # Murnaghan Equation of State Information
# #           parameters = (E0, K0, K0', V0)

#   Diamond structure calculations
initial_parameters_diamond = (total_energies_strain_diamond[mid(total_energies_strain_diamond)], 1e11, 3.5, volumes_sim_diamond[mid(volumes_sim_diamond)])
fit_parameters_diamond = sp.fmin(square_differences, initial_parameters_diamond, args=(volumes_sim_diamond, total_energies_strain_diamond, murnaghan), maxiter=100000)
volumes_diamond = np.linspace(volumes_sim_diamond[0], volumes_sim_diamond[-1], num=number_of_volume_points)
fit_total_energies_strain_diamond = murnaghan(fit_parameters_diamond, volumes_diamond)

#     Beta-Sn structure calculations
initial_parameters_betasn = (total_energies_strain_betasn[mid(total_energies_strain_betasn)], 1e11, 3.5, volumes_sim_betasn[mid(volumes_sim_betasn)])
fit_parameters_betasn = sp.fmin(square_differences, initial_parameters_betasn, args=(volumes_sim_betasn, total_energies_strain_betasn, murnaghan), maxiter=100000)
volumes_betasn = np.linspace(volumes_sim_betasn[0], volumes_sim_betasn[-1], num=number_of_volume_points)
fit_total_energies_strain_betasn = murnaghan(fit_parameters_betasn, volumes_betasn)

vol_mash = np.concatenate((volumes_sim_diamond, volumes_sim_betasn))
volume = np.linspace(min(vol_mash), max(vol_mash), num=number_of_enthalpy_samples)
pressure = np.linspace(-pressure_extreme, pressure_extreme, number_of_enthalpy_samples)
enthalpy_diamond = enthal_murn(fit_parameters_diamond, pressure)
enthalpy_beta = enthal_murn(fit_parameters_betasn, pressure)
matches = np.zeros(number_of_enthalpy_samples)
index = np.arange(0, number_of_enthalpy_samples, 1)

for x in index:
    add = np.array([enthalpy_diamond[x]/enthalpy_beta[x]])
    matches[x] = add

enthalpy_diamond = enthalpy_diamond[np.logical_not(np.isnan(matches))]
enthalpy_beta = enthalpy_beta[np.logical_not(np.isnan(matches))]
pressure = pressure[np.logical_not(np.isnan(matches))]
matches = matches[np.logical_not(np.isnan(matches))]
near = find_nearest(matches, 1)


tpressure = pressure[near[1]]
tvol_diamond = transition_volume(fit_parameters_diamond, tpressure)
tvol_beta = transition_volume(fit_parameters_betasn, tpressure)

#Parameters modifyied for wiki
vol_0_diamond = round(1e30*fit_parameters_diamond[3], 2)
k_0_diamond = round(1e-9*fit_parameters_diamond[1], 1)
k_0_prime_diamond = round(fit_parameters_diamond[2], 2)
vol_t_diamond = round(1e30*tvol_diamond, 2)
vol_0_betasn = round(1e30*fit_parameters_betasn[3], 2)
# celldm_ratio_betasn = round(celldm_3_betasn, 3)
k_0_betasn = round(1e-9*fit_parameters_betasn[1], 1)
k_0_prime_betasn = round(fit_parameters_betasn[2], 2)
vol_t_betasn = round(1e30*tvol_beta,2)
t_pres = round(1e-9*tpressure, 2)

# print(fit_parameters_diamond[0]+46.39167521*joules_per_rydberg)
print('diamond')
print(vol_0_diamond)
print(vol_t_diamond)
print(k_0_diamond)
print(k_0_prime_diamond)
print(t_pres)
print('betasn')
print(vol_0_betasn)
print(vol_t_betasn)
print(k_0_betasn)
print(k_0_prime_betasn)
print(t_pres)

# plt.plot(volumes_betasn,fit_total_energies_strain_betasn)
# plt.plot(volumes_diamond, fit_total_energies_strain_diamond)
# plt.show()