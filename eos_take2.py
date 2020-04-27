import scipy.optimize as sp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Use TeX fonts
# mpl.rcParams['text.usetex'] = Trueeos_take2.py
# mpl.rcParams['font.sans-serif'] = "cmr10"

# Parameters
number_of_volume_points = 100
number_of_enthalpy_samples = 10000
exchange_correlation = 'SCAN'
chemical_formula = 'Si'
structure_names = ['$Fd\overline{3}m$', '$I4_{1}/amd$']  # Hermann-Mauguin notation for diamond & beta-Sn structures
n_atom_diamond = 8
n_atom_betasn = 8

# Conversion factors
cubic_meters_per_cubic_angstrom = 1e-30
joules_per_Rydberg = 2.1798741e-18
# print(exchange_correlation)
# Data
if exchange_correlation == 'PBE':
    figure_file_name = 'Si.PBE_1ML.EoS.png'
    volumes_sim_BetaSn = cubic_meters_per_cubic_angstrom*np.array([9.65158069550,10.29501940960,10.93845812370,11.58189683782,12.86877426200,14.15565168618,15.44252911440,16.72940654262])
    total_energies_strain_BetaSn = (joules_per_Rydberg/n_atom_betasn)*np.array([-372.66327325,-373.06805755,-373.30465661,-373.43427048,-373.50849372,-373.47259800,-373.40956222,-373.34922248])
    volumes_sim_diamond = cubic_meters_per_cubic_angstrom*np.array([14.73811721082,16.84356252521,18.94900784461,21.05445315901,23.15989847344,25.26534379280,27.37078910724])
    total_energies_strain_diamond = (joules_per_Rydberg/n_atom_betasn)*np.array([-373.29868028,-373.41432946,-373.48607830,-373.51022345,-373.48973472,-373.44785756,-373.38548468])
    ########older that 4/25/2020 data
    # volumes_sim_diamond = cubic_meters_per_cubic_angstrom*np.array([14.73811721082,16.84356252521,18.94900784461,20.63336409813,20.84390862609,21.05445315901,21.05445315901,21.26499769197,21.47554221992,23.15989847344,25.26534379280,27.37078910724])
    # total_energies_strain_diamond = (joules_per_Rydberg/n_atom_diamond)*np.array([-373.29868028,-373.41432946,-373.48607830,-373.50907404,-373.50991929,-373.51022345,-373.51022345,-373.50988601,-373.50908026,-373.48973472,-373.44785756,-373.38548468])
    # #   Original data
    # # volumes_sim_BetaSn = cubic_meters_per_cubic_angstrom*np.array([18.37914863026,21.00474128854,23.63033394677,25.73080807936,25.99336734596,26.25592661255,26.25592661255,26.51848587909,26.78104514568,28.88151927827,31.50711193656,34.13270459479])
    # # total_energies_strain_BetaSn = (joules_per_Rydberg/n_atom_betasn)*np.array([-186.55130659,-186.54211094,-186.54976560,-186.55931303,-186.55961402,-186.55971389,-186.55971389,-186.55961507,-186.55931838,-186.55006955,-186.52439387,-186.47946622])
    # #
    # #   First data point removed
    # # volumes_sim_BetaSn = cubic_meters_per_cubic_angstrom*np.array([21.00474128854,23.63033394677,25.73080807936,25.99336734596,26.25592661255,26.25592661255,26.51848587909,26.78104514568,28.88151927827,31.50711193656,34.13270459479])
    # # total_energies_strain_BetaSn = (joules_per_Rydberg/n_atom_betasn)*np.array([-186.54211094,-186.54976560,-186.55931303,-186.55961402,-186.55971389,-186.55971389,-186.55961507,-186.55931838,-186.55006955,-186.52439387,-186.47946622])
    # #
    # #   First two data points removed
    # volumes_sim_BetaSn = cubic_meters_per_cubic_angstrom*np.array([23.63033394677,25.73080807936,25.99336734596,26.25592661255,26.25592661255,26.51848587909,26.78104514568,28.88151927827,31.50711193656,34.13270459479])
    # total_energies_strain_BetaSn = (joules_per_Rydberg/n_atom_betasn)*np.array([-186.54976560,-186.55931303,-186.55961402,-186.55971389,-186.55971389,-186.55961507,-186.55931838,-186.55006955,-186.52439387,-186.47946622])
if exchange_correlation == 'PZ':
    figure_file_name = 'Si.PZ_1ML.EoS.png'
    volumes_sim_BetaSn = cubic_meters_per_cubic_angstrom*np.array([8.69553212657,9.93775100516,11.17996987587,12.42218875449,13.66440763310,14.90662650382,16.14884538240])
    total_energies_strain_BetaSn = (joules_per_Rydberg/n_atom_betasn)*np.array([-89.22482893,-90.30796765,-90.68395156,-90.76086191,-90.72389896,-90.65947349,-90.59691208])
    volumes_sim_diamond = cubic_meters_per_cubic_angstrom*np.array([14.08392451545,16.09591373539,18.10790294570,20.11989216564,22.13188138559,24.14387059586,26.15585981581])
    total_energies_strain_diamond = (joules_per_Rydberg/n_atom_betasn)*np.array([-90.51829730,-90.62944720,-90.70666769,-90.73095209,-90.70955754,-90.66993654,-90.60806563])
    ###########older than 4/25/2020
    # volumes_sim_diamond = cubic_meters_per_cubic_angstrom*np.array([14.08392451545,16.09591373539,18.10790294570,19.71749432455,19.91869324267,20.11989216564,20.11989216564,20.32109108858,20.52229000673,22.13188138559,24.14387059586,26.15585981581])
    # total_energies_strain_diamond = (joules_per_Rydberg/n_atom_diamond)*np.array([-90.51829730,-90.62944720,-90.70666769,-90.72996962,-90.73070834,-90.73095209,-90.73095209,-90.73070994,-90.72997869,-90.70955754,-90.66993654,-90.60806563])
    # #   Original data
    # # volumes_sim_BetaSn = cubic_meters_per_cubic_angstrom*np.array([11.70887675621,13.38157343251,15.05427010886,16.39242745286,16.55969712631,16.72696679248,16.72696679248,16.89423645866,17.06150613210,18.39966347610,20.07236015245,21.74505682875])
    # # total_energies_strain_BetaSn = (joules_per_Rydberg/n_atom_betasn)*np.array([-45.21547853,-45.19559075,-45.21462844,-45.21920332,-45.21922258,-45.21913143,-45.21913143,-45.21893224,-45.21862437,-45.21185813,-45.19217933,-45.16288506])
    # #
    # #   First data point removed
    # volumes_sim_BetaSn = cubic_meters_per_cubic_angstrom*np.array([13.38157343251,15.05427010886,16.39242745286,16.55969712631,16.72696679248,16.72696679248,16.89423645866,17.06150613210,18.39966347610,20.07236015245,21.74505682875])
    # total_energies_strain_BetaSn = (joules_per_Rydberg/n_atom_betasn)*np.array([-45.19559075,-45.21462844,-45.21920332,-45.21922258,-45.21913143,-45.21913143,-45.21893224,-45.21862437,-45.21185813,-45.19217933,-45.16288506])

if exchange_correlation == 'SCAN':
    figure_file_name = 'Si.SCAN_1ML.EoS.png'
    volumes_sim_BetaSn = cubic_meters_per_cubic_angstrom*np.array([10.21619778257,11.67565460811,13.13511143364,14.59456825921,16.05402508474,17.51348191028,18.97293873584])
    total_energies_strain_BetaSn = (joules_per_Rydberg/n_atom_betasn)*np.array([-61.24625462,-62.21356102,-62.54424616,-62.61104424,-62.57780653,-62.51573929,-62.45101265])
    volumes_sim_diamond = cubic_meters_per_cubic_angstrom*np.array([13.94338719066,15.93529964919,17.92721209824,19.91912455677,21.91103701526,23.90294946430,25.89486192283])
    total_energies_strain_diamond = (joules_per_Rydberg/n_atom_diamond)*np.array([-62.34036437,-62.49221543,-62.53738742,-62.53699016,-62.50937811,-62.43660083,-62.35068451])
    ############older tha 4/25/2020
    # volumes_sim_diamond = cubic_meters_per_cubic_angstrom*np.array([15.93529964919,17.92721209824,19.52074206884,19.71993331044,19.91912455677,19.91912455677,20.11831580309,20.31750704466,21.91103701526,23.90294946430])
    # total_energies_strain_diamond = (joules_per_Rydberg/n_atom_diamond)*np.array([-62.34036437,-62.49221543,-62.53738742,-62.53589659,-62.53676157,-62.53699016,-62.53684557,-62.53579955,-62.50937811,-62.43660083])
    # #   Original data
    # # volumes_sim_BetaSn = cubic_meters_per_cubic_angstrom*np.array([17.67855882259,20.20406722585,22.72957562906,24.74998235166,25.00253319197,25.25508403227,25.25508403227,25.50763487263,25.76018571294,27.78059243554,30.30610083875,32.83160924196])
    # # total_energies_strain_BetaSn = (joules_per_Rydberg/n_atom_betasn)*np.array([-31.11700645,-31.06535235,-31.08317680,-31.09456583,-31.09493016,-31.09515189,-31.09515189,-31.09511682,-31.09486093,-31.08511596,-31.05350619,-31.00657811])
    # #   First data point removed
    # volumes_sim_BetaSn = cubic_meters_per_cubic_angstrom*np.array([20.20406722585,22.72957562906,24.74998235166,25.00253319197,25.25508403227,25.25508403227,25.50763487263,25.76018571294,27.78059243554,30.30610083875,32.83160924196])
    # total_energies_strain_BetaSn = (joules_per_Rydberg/n_atom_betasn)*np.array([-31.06535235,-31.08317680,-31.09456583,-31.09493016,-31.09515189,-31.09515189,-31.09511682,-31.09486093,-31.08511596,-31.05350619,-31.00657811])

def mid(x):
    return int((len(x)-1)/2)


def square_differences(p, x, y, f):
    return np.power(f(p, x) - y, 2).sum()


def murnaghan(p, v):
    kk = p[2]-1.0
    return (p[0] + (p[1]*p[3]*(((1.0/(p[2]*kk))*np.power((v/p[3]), (-kk)))+(v/(p[2]*p[3]))-(1.0/kk))))


# Murnaghan Equation of State Information
#           parameters = (E0, K0, K0', V0)

#   Diamond structure calculations
initial_parameters_diamond = (total_energies_strain_diamond[mid(total_energies_strain_diamond)], 1e11, 3.5, volumes_sim_diamond[mid(volumes_sim_diamond)])
fit_parameters_diamond = sp.fmin(square_differences, initial_parameters_diamond, args=(volumes_sim_diamond, total_energies_strain_diamond, murnaghan), maxiter=100000)
volumes_diamond = np.linspace(volumes_sim_diamond[0], volumes_sim_diamond[-1], num=number_of_volume_points)
fit_total_energies_strain_diamond = murnaghan(fit_parameters_diamond, volumes_diamond)

#     Beta-Sn structure calculations
volumes_BetaSn = np.linspace(volumes_sim_BetaSn[0], volumes_sim_BetaSn[-1], num=number_of_volume_points)
#
        #cell_dofree ='z' values
initial_parameters_shape_BetaSn = (total_energies_strain_BetaSn[mid(total_energies_strain_BetaSn)], 1e11, 3.5, volumes_sim_BetaSn[mid(volumes_sim_BetaSn)])
fit_parameters_shape_BetaSn = sp.fmin(square_differences, initial_parameters_shape_BetaSn, args=(volumes_sim_BetaSn, total_energies_strain_BetaSn, murnaghan), maxiter=100000)
fit_total_energies_strain_BetaSn = murnaghan(fit_parameters_shape_BetaSn, volumes_BetaSn)
# print(fit_parameters_z_BetaSn)


# Phase transition properties
#       parameters = (E0, K0, K0', V0)
# Murnaghan equation-of-state enthalpy
def enthal_murn(a, p):
    e0 = a[0]
    k0 = a[1]
    k0prime = a[2]
    v0 = a[3]
    if (1+p*(k0prime/k0)).any() > 0 and k0 != 0. and k0prime != 0.:
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
pressure = np.linspace(-50*11.7e9, 50*11.7e9, number_of_enthalpy_samples)
diamond = enthal_murn(fit_parameters_diamond, pressure)
beta = enthal_murn(fit_parameters_shape_BetaSn, pressure)
matches = np.zeros(number_of_enthalpy_samples)
index = np.arange(0,number_of_enthalpy_samples,1)



#find indexes were plots intersect
enthal_inter_idx=np.argwhere(np.diff(np.nan_to_num(np.sign(beta-diamond))))
n=0
for i in enthal_inter_idx:
    if np.isnan(beta[i]):
        enthal_inter_idx=np.delete(enthal_inter_idx, n)

#chose either positive or negative transition pressure
diamond_energy_min=np.amin(fit_total_energies_strain_diamond)
betasn_energy_min=np.amin(fit_parameters_shape_BetaSn)
energy_difference=-diamond_energy_min+betasn_energy_min
# for i in enthal_inter_idx:
#     signofp=int(np.sign(pressure[i]))
#     # print(pressure[i])
#     if signofp==int(np.sign(energy_difference)):
#         tpressure=pressure[i]

tpressure=np.amin(pressure[enthal_inter_idx])
for x in index:
    add = np.array([diamond[x]/beta[x]])
    matches[x] = add

diamond = diamond[np.logical_not(np.isnan(matches))]
beta = beta[np.logical_not(np.isnan(matches))]
pressure = pressure[np.logical_not(np.isnan(matches))]
matches = matches[np.logical_not(np.isnan(matches))]
near = find_nearest(matches, 1)

# tpressure = pressure[near[1]]
tvol_diamond = transition_volume(fit_parameters_diamond, tpressure)
tvol_beta = transition_volume(fit_parameters_shape_BetaSn, tpressure)

#Parameters modifyied for wiki
vol_0_diamond = round(1e30*fit_parameters_diamond[3], 2)
k_0_diamond = round(1e-9*fit_parameters_diamond[1], 1)
k_0_prime_diamond = round(fit_parameters_diamond[2], 2)
vol_t_diamond = round(1e30*tvol_diamond, 2)
vol_0_betasn = round(1e30*fit_parameters_shape_BetaSn[3], 2)
# celldm_ratio_betasn = round(celldm_3_BetaSn, 3)
k_0_betasn = round(1e-9*fit_parameters_shape_BetaSn[1], 1)
k_0_prime_betasn = round(fit_parameters_shape_BetaSn[2], 2)
vol_t_betasn = round(1e30*tvol_beta,2)
t_pres = round(1e-9*tpressure, 2)

# print(murnaghan(fit_parameters_diamond,tvol_diamond))
print(tpressure*1e-9)
print(tvol_diamond/cubic_meters_per_cubic_angstrom)
print(tvol_beta/cubic_meters_per_cubic_angstrom)

# # Two-phase plot
# fig9 = plt.figure()
# plt.plot(volume/cubic_meters_per_cubic_angstrom, (murnaghan(fit_parameters_shape_BetaSn, tvol_beta)-(volume-tvol_beta)*tpressure)*6.242e18,color='y')
# plt.plot(volumes_diamond/cubic_meters_per_cubic_angstrom, fit_total_energies_strain_diamond*6.242e18,color='slateblue')
# plt.scatter(volumes_sim_diamond/cubic_meters_per_cubic_angstrom, total_energies_strain_diamond*6.242e18,label=r'Diamond',color='slateblue',marker='s')
# plt.plot(volumes_BetaSn/cubic_meters_per_cubic_angstrom, fit_total_energies_strain_BetaSn*6.242e18,color='firebrick')
# plt.scatter(volumes_sim_BetaSn/cubic_meters_per_cubic_angstrom, total_energies_strain_BetaSn*6.242e18, label=r'BetaSn',color='firebrick',marker='v')
# plt.yscale(r'linear')
# plt.title(r'Energy Equation of State for ' + structure_names[0] + ' and ' + structure_names[1] + ' ' +
#                chemical_formula)
# plt.ylabel(r'Total Energy (eV/atom)')
# plt.xlabel(r'Volume (\r{A}$^3$/atom)')
# plt.legend()
#
# # ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax))
# # plt.plot([vol_0_diamond,vol_0_diamond],[-180,murnaghan(fit_parameters_diamond,vol_0_diamond)])
# # plt.text(vol_0_diamond, 0,'V_0',verticalalignment='center')
#
#
# plt.show()
# print(t_pres)
# print(tvol_diamond)
# print(tvol_beta)
# print(vol_0_betasn)