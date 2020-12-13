import scipy.optimize as sp
import matplotlib as mpl

# Parameters
number_of_volume_points = 100
number_of_enthalpy_samples = 100000
pressure_extreme = 50*11.7e9
exchange_correlation = 'SCAN'
if exchange_correlation == 'PZ':
    from QEData_pz import *
elif exchange_correlation == 'PBE':
    from QEData import *
    figure_file_name = 'Si.PBE.EoS.png'
elif exchange_correlation == 'SCAN':
    from QEData_scan import *
    # from QEData_scan_paper import *

# Conversion factors
# from unit_conversions import cubic_meters_per_cubic_angstrom, joules_per_rydberg

figure_file_name = 'Si.' + exchange_correlation + '.EoS.png'

# Use TeX fonts
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.sans-serif'] = "cmr10"

def mid(x):
    return int((len(x)-1)/2)


def square_differences(p, x, y, f):
    return np.power(f(p, x) - y, 2).sum()


def murnaghan(p, v):
    kk = p[2]-1.0
    return p[0] + (p[1] * p[3] * (((1.0 / (p[2] * kk)) * np.power((v / p[3]), (-kk))) +
                                  (v / (p[2] * p[3])) - (1.0 / kk)))


# Murnaghan Equation of State Information
#           parameters = (E0, K0, K0', V0)

#   Diamond structure calculations
volumes_sim_diamond = np.power(lattice_parameters_diamond, 3) / n_atom_diamond
initial_parameters_diamond = (total_energies_strain_diamond[mid(total_energies_strain_diamond)], 1e11, 3.5, volumes_sim_diamond[mid(volumes_sim_diamond)])
fit_parameters_diamond = sp.fmin(square_differences, initial_parameters_diamond, args=(volumes_sim_diamond, total_energies_strain_diamond, murnaghan), maxiter=100000)
volumes_diamond = np.linspace(volumes_sim_diamond[0], volumes_sim_diamond[-1], num=number_of_volume_points)
fit_total_energies_strain_diamond = murnaghan(fit_parameters_diamond, volumes_diamond)

#   Beta-Sn structure calculations
volumes_sim_BetaSn = (celldm_3_BetaSn / n_atom_betasn)*np.power(lattice_parameters_BetaSn, 3)
volumes_BetaSn = np.linspace(volumes_sim_BetaSn[0], volumes_sim_BetaSn[-1], num=number_of_volume_points)

#       cell_dofree ='all' values
initial_parameters_all_BetaSn = (total_energies_strain_all_BetaSn[mid(total_energies_strain_all_BetaSn)], 1e11, 3.5, volumes_sim_BetaSn[mid(volumes_sim_BetaSn)])
fit_parameters_all_BetaSn = sp.fmin(square_differences, initial_parameters_all_BetaSn, args=(volumes_sim_BetaSn, total_energies_strain_all_BetaSn, murnaghan), maxiter=100000)
fit_total_energies_strain_all_BetaSn = murnaghan(fit_parameters_all_BetaSn, volumes_BetaSn)
# print(fit_parameters_all_BetaSn)

#       cell_dofree ='z' values
initial_parameters_shape_BetaSn = (total_energies_strain_shape_BetaSn[mid(total_energies_strain_shape_BetaSn)], 1e11, 3.5, volumes_sim_BetaSn[mid(volumes_sim_BetaSn)])
fit_parameters_shape_BetaSn = sp.fmin(square_differences, initial_parameters_shape_BetaSn, args=(volumes_sim_BetaSn, total_energies_strain_shape_BetaSn, murnaghan), maxiter=100000)
fit_total_energies_strain_shape_BetaSn = murnaghan(fit_parameters_shape_BetaSn, volumes_BetaSn)
# print(fit_parameters_z_BetaSn)


# Phase transition properties
#       parameters = (E0, K0, K0', V0)
# Murnaghan equation-of-state enthalpy
def enthal_murn(a, p):
    e0 = a[0]
    k0 = a[1]
    k0prime = a[2]
    v0 = a[3]
    if (1+p*(k0prime/k0)).any() > 0 and k0 != 0 and k0prime != 0:
        v = v0*np.power((1 + p * (k0prime / k0)), (-1.0 / k0prime))
    enthalpy = e0
    enthalpy += (k0 * v0 /(k0prime * (k0prime - 1.0))) * np.power(v/v0, (1 - k0prime))
    enthalpy += (k0 * v0 / k0prime) * (v/v0)
    enthalpy += -(k0 * v0 /(k0prime - 1))
    enthalpy += p*v
    return enthalpy


# Find nearest neighbor atoms
def find_nearest(array, value):
    array = np.asarray(array)
    idx = np.argmin((np.abs(array - value)))
    return [array[idx], idx]


# Transition volume
def transition_volume(a, p):
    e0 = a[0]
    k0 = a[1]
    k0prime = a[2]
    v0 = a[3]
    v = v0 * np.power((1 + p * (k0prime / k0)), (-1.0 / k0prime))
    return v


vol_mash = np.concatenate((volumes_sim_diamond, volumes_sim_BetaSn))
volume = np.linspace(min(vol_mash), max(vol_mash), num=number_of_enthalpy_samples)
pressure = np.linspace(-pressure_extreme, pressure_extreme, number_of_enthalpy_samples)
diamond = enthal_murn(fit_parameters_diamond, pressure)
beta = enthal_murn(fit_parameters_shape_BetaSn, pressure)
matches = np.zeros(number_of_enthalpy_samples)
index = np.arange(0, number_of_enthalpy_samples, 1)

for x in index:
    add = np.array([diamond[x]/beta[x]])
    matches[x] = add

diamond = diamond[np.logical_not(np.isnan(matches))]
beta = beta[np.logical_not(np.isnan(matches))]
pressure = pressure[np.logical_not(np.isnan(matches))]
matches = matches[np.logical_not(np.isnan(matches))]
near = find_nearest(matches, 1)


tpressure = pressure[near[1]]
tvol_diamond = transition_volume(fit_parameters_diamond, tpressure)
tvol_beta = transition_volume(fit_parameters_shape_BetaSn, tpressure)

#Parameters modifyied for wiki
vol_0_diamond = round(1e30*fit_parameters_diamond[3], 2)
k_0_diamond = round(1e-9*fit_parameters_diamond[1], 1)
k_0_prime_diamond = round(fit_parameters_diamond[2], 2)
vol_t_diamond = round(1e30*tvol_diamond, 2)
vol_0_betasn = round(1e30*fit_parameters_shape_BetaSn[3], 2)
celldm_ratio_betasn = round(celldm_3_BetaSn, 3)
k_0_betasn = round(1e-9*fit_parameters_shape_BetaSn[1], 1)
k_0_prime_betasn = round(fit_parameters_shape_BetaSn[2], 2)
vol_t_betasn = round(1e30*tvol_beta,2)
t_pres = round(1e-9*tpressure, 2)

# print(fit_parameters_diamond[0]+46.39167521*joules_per_rydberg)
print(vol_0_diamond)
print(vol_t_diamond)
print(t_pres)
print(vol_0_betasn)
print(vol_t_betasn)
print(k_0_betasn)
print(k_0_diamond)
print(fit_parameters_diamond[0]/1.6e-19)
print(fit_parameters_shape_BetaSn[0]/1.6e-19)