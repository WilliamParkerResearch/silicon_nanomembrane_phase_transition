from ReferenceFiles.equations_of_state import *
from ReferenceFiles.FunctionDefinitions import *
from ReferenceFiles.enthalpy import *
import matplotlib.pyplot as plt
import scipy.optimize as sp
from scipy import interpolate
import numpy as np

number_of_volume_points = 1000000
n_pressure_values = 1000000
exchange_correlation = 'PBE'
N_ML = '2'

directoryofdata = 'DataFolder' + '.' + exchange_correlation + '.' + 'eosm' + '.' + 'Data_' + N_ML + 'L'
exec('from {} import *'.format(directoryofdata))

eos_type = 'murnaghan'

if eos_type == 'birch-murnaghan':
    def eos(p, v):
        return birch_murnaghan(p, v)
elif eos_type == 'murnaghan':
    def eos(p, v):
        return murnaghan(p, v)
elif eos_type == 'vinet':
    def eos(p, v):
        return vinet(p, v)


def intersection_idx(equation):  # prints lowerbound index in which the value is found within it and the index above it
    i = np.argwhere(np.diff(np.sign(equation))).flatten()
    return i

def extend_array_extremas(array,percentage=0):
    min_value = np.amin(array)
    max_value = np.amax(array)
    value_diff = (percentage)*np.absolute(max_value-min_value)
    return value_diff

#   Diamond structure calculations
initial_parameters_diamond = (total_energies_strain_diamond[mid(total_energies_strain_diamond)], 1e11, 3.5,volumes_sim_diamond[mid(volumes_sim_diamond)])
fit_parameters_diamond = sp.fmin(square_differences, initial_parameters_diamond, args=(volumes_sim_diamond, total_energies_strain_diamond, eos), maxiter=100000)
volume_diamond_extension = extend_array_extremas(volumes_sim_diamond,0.5)
volumes_diamond = np.linspace(np.amin(volumes_sim_diamond)-volume_diamond_extension, np.amax(volumes_sim_diamond)+volume_diamond_extension, num=number_of_volume_points)
pressures_diamond = -pressure_from_energy_equation_of_state(fit_parameters_diamond,volumes_diamond,eos=eos_type)
enthalpy_diamond = enthalpy_from_volume(fit_parameters_diamond, volumes_diamond, eos=eos_type)
fit_total_energies_strain_diamond = eos(fit_parameters_diamond, volumes_diamond)
#     Beta-Sn structure calculations
initial_parameters_shape_betasn = (total_energies_strain_betasn[mid(total_energies_strain_betasn)], 1e11, 3.5, volumes_sim_betasn[mid(volumes_sim_betasn)])
fit_parameters_betasn = sp.fmin(square_differences, initial_parameters_shape_betasn, args=(volumes_sim_betasn, total_energies_strain_betasn, eos), maxiter=100000)
volume_betasn_extension = extend_array_extremas(volumes_sim_betasn,0.5)
volumes_betasn = np.linspace(np.amin(volumes_sim_betasn)-volume_betasn_extension, np.amax(volumes_sim_betasn)+volume_betasn_extension, num=number_of_volume_points)
pressures_betasn = -pressure_from_energy_equation_of_state(fit_parameters_betasn,volumes_betasn,eos=eos_type)
enthalpy_betasn = enthalpy_from_volume(fit_parameters_betasn, volumes_betasn, eos=eos_type)
fit_total_energies_strain_betasn = eos(fit_parameters_betasn, volumes_betasn)
#           parameters = (E0, K0, K0', V0)


# plt.plot(pressures_betasn, enthalpy_betasn*1e18,color='r')
# plt.plot(pressures_diamond, enthalpy_diamond*1e18,color='b')
# plt.show()


idx_d = np.argsort(pressures_diamond)
pressures_diamond = pressures_diamond[idx_d]
enthalpy_diamond = enthalpy_diamond[idx_d]

idx_b = np.argsort(pressures_betasn)
pressures_betasn = pressures_betasn[idx_b]
enthalpy_betasn = enthalpy_betasn[idx_b]

pressure_min = np.amax(np.array([pressures_diamond[0],pressures_betasn[0]]))
pressure_max = np.amin(np.array([pressures_diamond[-1],pressures_betasn[-1]]))
pressures = np.linspace(pressure_min,pressure_max,n_pressure_values)

enthalpy_diamond_func = interpolate.interp1d(pressures_diamond,enthalpy_diamond)
enthalpy_diamond_fit = enthalpy_diamond_func(pressures)
enthalpy_betasn_func = interpolate.interp1d(pressures_betasn,enthalpy_betasn)
enthalpy_betasn_fit = enthalpy_betasn_func(pressures)

enthalpy_diff = enthalpy_diamond_fit-enthalpy_betasn_fit

p_idx = intersection_idx(enthalpy_diff)
t_pressure = pressures[int(p_idx)]

t_vol_d_idx = intersection_idx(pressures_diamond-t_pressure)
t_volume_diamond = volumes_diamond[t_vol_d_idx]

t_vol_b_idx = intersection_idx(pressures_betasn-t_pressure)
t_volume_betasn = volumes_betasn[t_vol_b_idx]


print(t_pressure*1e-9,'GPa')
print(t_volume_diamond*1e30,'Angstroms^3')
print(t_volume_betasn*1e30,'Angstroms^3')
#
# plt.plot(pressures, enthalpy_diamond_fit*1e18)
# plt.plot(pressures, enthalpy_betasn_fit*1e18)
# plt.plot(pressures_betasn, enthalpy_betasn*1e18)
# plt.plot(pressures_diamond, enthalpy_diamond*1e18)
# plt.xlim([pressure_min,pressure_max])
# plt.ylim([1e18*np.amin(np.append(enthalpy_diamond_fit,enthalpy_betasn_fit)),1e18*np.amax(np.append(enthalpy_diamond_fit,enthalpy_betasn_fit))])
# plt.show()
#
# plt.plot(volumes_diamond,pressures_diamond)
# plt.hlines(y=t_pressure,xmin=volumes_diamond[0],xmax=volumes_diamond[-1])
# plt.show()
#
# plt.plot(volumes_betasn,pressures_betasn)
# plt.hlines(y=t_pressure,xmin=volumes_betasn[0],xmax=volumes_betasn[-1])
# plt.show()
#
# all_volumes = np.append(volumes_betasn, volumes_diamond)
# volumes_range = np.linspace(np.amin(all_volumes),np.amax(all_volumes[-1]),number_of_volume_points)
# tangent = t_pressure * (volumes_range - t_volume_betasn)+eos(fit_parameters_betasn, t_volume_betasn)
#
# plt.plot(volumes_betasn, fit_total_energies_strain_betasn)
# plt.plot(volumes_diamond, fit_total_energies_strain_diamond)
# plt.plot(volumes_range, tangent)
# plt.show()