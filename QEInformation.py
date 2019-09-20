# Things to do:
#   Moving the individual graph properties (xlabels, ylabels, etc.) to arrays so we can easily modify their formatting
#   Fix the scale of the y-axes in certain plots,
#   Changing the conversion factors as variable names: ? = 6.242e18
#   Convert total energies in convergence plots to difference in total energy from converged total energy
#   Put Hermann-Mauguin notation in plot legends
#
import scipy.optimize as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
from QEData import *


# Use TeX fonts
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.sans-serif'] = "cmr10"

# Simulation details
chemical_formula = 'Si'
structure_names = ['$Fd\overline{3}m$', '$I4_{1}/amd$']  # Hermann-Mauguin notation for diamond & beta-Sn structures
number_of_atoms = [8, 4]
convergence_threshold = 7.35e-5     # in Ry = 1 meV

# Graph properties
plot_titles = [r'Cutoff Energy Convergence of Total Energy for ' + structure_names[0] + ' ' + chemical_formula,
               r'Cutoff Energy Convergence of Total Energy for ' + structure_names[1] + ' ' + chemical_formula,
               r'Cutoff Energy Convergence of ' + structure_names[0] + ' and ' + structure_names[1] + ' ' +
               chemical_formula + ' Total Energy Difference',
               r'$\vec{k}$-point Mesh Convergence of Total Energy for ' + structure_names[0] + ' ' + chemical_formula,
               r'$\vec{k}$-point Mesh Convergence of Total Energy for ' + structure_names[1] + ' ' + chemical_formula,
               r'$\vec{k}$-point Mesh Convergence of ' + structure_names[0] + ' and ' + structure_names[1] + ' ' +
               chemical_formula + ' Total Energy Difference',
               r'Energy Equation of State for ' + structure_names[0] + ' ' + chemical_formula,
               r'Energy Equation of State for ' + structure_names[1] + ' ' + chemical_formula,
               r'Energy Equation of State for ' + structure_names[0] + ' and ' + structure_names[1] + ' ' +
               chemical_formula
               ]

plot_xlabels = [r"Cutoff Energy (eV)",
                r"Cutoff Energy (eV)",
                r"Cutoff Energy (eV)",
                r"Number of $\vec{k}$ points in each direction",
                r"Number of $\vec{k}$ points in each direction",
                r"Number of $\vec{k}$ points in each direction",
                r'Volume (m$^3$/atom)'
                ]

plot_ylabels = []

# Unit conversions
cubic_meters_per_cubic_angstrom = 1e-30
pascals_per_gigapascal = 1e9


# Function definitions
def round(n, decimals=0):
    multiplier = 10 ** decimals
    return int(n * multiplier) / multiplier


def murnaghan(p, v):
    kk = p[2]-1.0
    return (p[0] + (p[1]*p[3]*(((1.0/(p[2]*kk))*np.power((v/p[3]), (-kk)))+(v/(p[2]*p[3]))-(1.0/kk))))


def square_differences(p, x, y, f):
    return np.power(f(p, x) - y, 2).sum()


def mid(x):
    return int((len(x)-1)/2)


def matchmaker(a, b):
    list_a = list(a)
    list_b = list(b)

    matches = np.intersect1d(a, b)
    index_a = []
    for x in matches:
        index_a.append(list_a.index(x))
    index_b = []
    for x in matches:
        index_b.append(list_b.index(x))
    return [index_a, index_b]


# Cutoff energy plots
#    Diamond structure
figure_index = 0
fig1 = plt.figure(figure_index+1)
plt.plot(cutoff_energies_diamond, total_energies_ecut_diamond, label='Total energy')
plt.scatter(cutoff_energies_diamond, total_energies_ecut_diamond)
plt.title(plot_titles[figure_index])
plt.xlabel(plot_xlabels[figure_index])
plt.ylabel(r"Total Energy (meV/atom)")
plt.axhline(total_energies_ecut_diamond[-1], color='r', label='convergence energy')
plt.axhline(total_energies_ecut_diamond[-1]-convergence_threshold, color='r', linestyle='--',
            label='convergence threshold')
plt.axhline(total_energies_ecut_diamond[-1]+convergence_threshold, color='r', linestyle='--')
plt.legend()
# plt.axis([0, ecut[-1],etot[-1]-10*convergence_threshold, etot[-1]+10*convergence_threshold])


#   Beta-Sn structure
figure_index += 1
fig2 = plt.figure(figure_index+1)
plt.plot(cutoff_energies_BetaSn, total_energies_ecut_BetaSn, label='Total energy')
plt.scatter(cutoff_energies_BetaSn, total_energies_ecut_BetaSn)
plt.title(plot_titles[figure_index])
plt.xlabel(plot_xlabels[figure_index])
plt.ylabel("Total Energy (meV/atom)")
plt.axhline(total_energies_ecut_BetaSn[-1], color='r', label='convergence energy')
plt.axhline(total_energies_ecut_BetaSn[-1]-convergence_threshold, color='r', linestyle='--',
            label='convergence threshold')
plt.axhline(total_energies_ecut_BetaSn[-1]+convergence_threshold, color='r', linestyle='--')
plt.legend()
# plt.axis([0, ecut[-1],etot[-1]-10*convergence_threshold, etot[-1]+10*convergence_threshold])

#   Difference between diamond and beta-Sn structures
indexes_ecut_diamond_BetaSn = matchmaker(cutoff_energies_diamond, cutoff_energies_BetaSn)
total_energies_ecut_diamond_matched = total_energies_kpoint_diamond[indexes_ecut_diamond_BetaSn[0]]
total_energies_ecut_BetaSn_matched = total_energies_kpoint_BetaSn[indexes_ecut_diamond_BetaSn[1]]
cutoff_energies_diamond_matched = cutoff_energies_diamond[indexes_ecut_diamond_BetaSn[0]]
cutoff_energies_BetaSn_matched = cutoff_energies_BetaSn[indexes_ecut_diamond_BetaSn[1]]

total_energies_ecut_diff = np.subtract(total_energies_ecut_BetaSn_matched, total_energies_ecut_diamond_matched)

figure_index += 1
fig3 = plt.figure(figure_index+1)
plt.plot(cutoff_energies_diamond_matched,total_energies_ecut_diff,label='energy difference')
plt.scatter(cutoff_energies_diamond_matched, total_energies_ecut_diff)
plt.title(plot_titles[figure_index])
plt.xlabel(plot_xlabels[figure_index])
plt.ylabel(r"Total Energy Difference (meV/atom)")


# k-point plots
#   Diamond structure
figure_index += 1
fig4 = plt.figure(figure_index+1)
plt.plot(kpoints_diamond,total_energies_kpoint_diamond, label='total energy')
plt.scatter(kpoints_diamond, total_energies_kpoint_diamond)
plt.title(plot_titles[figure_index])
plt.xlabel(plot_xlabels[figure_index])
plt.ylabel(r"Total Energy (meV/atom)")
plt.axhline(total_energies_kpoint_diamond[-1], color='r', label='convergence energy')
plt.axhline(total_energies_kpoint_diamond[-1]-convergence_threshold, color='r', linestyle='--', label='convergence threshold')
plt.axhline(total_energies_kpoint_diamond[-1]+convergence_threshold, color='r', linestyle='--')
plt.axis([0, kpoints_diamond[-1],total_energies_kpoint_diamond[-1]-10*convergence_threshold, total_energies_kpoint_diamond[-1]+10*convergence_threshold])
plt.legend()

#   beta-Sn structure
figure_index += 1
fig5 = plt.figure(figure_index+1)
plt.plot(kpoints_BetaSn,total_energies_kpoint_BetaSn, label='total energy')
plt.scatter(kpoints_BetaSn, total_energies_kpoint_BetaSn)
plt.title(plot_titles[figure_index])
plt.xlabel(plot_xlabels[figure_index])
plt.ylabel(r"Total Energy (meV/atom)")
plt.axhline(total_energies_kpoint_BetaSn[-1], color='r', label='convergence energy')
plt.axhline(total_energies_kpoint_BetaSn[-1]-convergence_threshold, color='r', linestyle='--', label='convergence threshold')
plt.axhline(total_energies_kpoint_BetaSn[-1]+convergence_threshold, color='r', linestyle='--')
plt.axis([0, kpoints_BetaSn[-1],total_energies_kpoint_BetaSn[-1]-10*convergence_threshold, total_energies_kpoint_BetaSn[-1]+10*convergence_threshold])
plt.legend()

#   Difference Between diamond and beta-Sn structures
indexes_kpoint_diamond_BetaSn = matchmaker(kpoints_diamond, kpoints_BetaSn)
total_energies_kpoint_diamond_matched = total_energies_kpoint_diamond[indexes_kpoint_diamond_BetaSn[0]]
total_energies_kpoint_BetaSn_matched = total_energies_kpoint_BetaSn[indexes_kpoint_diamond_BetaSn[1]]
kpoints_diamond_matched = kpoints_diamond[indexes_kpoint_diamond_BetaSn[0]]
kpoints_BetaSn_matched = kpoints_BetaSn[indexes_kpoint_diamond_BetaSn[1]]

total_energies_kpoint_diff = np.subtract(total_energies_kpoint_BetaSn_matched, total_energies_kpoint_diamond_matched)

figure_index += 1
fig6 = plt.figure(figure_index+1)
plt.plot(kpoints_diamond_matched,total_energies_kpoint_diff)
plt.scatter(kpoints_diamond_matched, total_energies_kpoint_diff)
plt.title(plot_titles[figure_index])
plt.xlabel(plot_xlabels[figure_index])
plt.ylabel(r"Total Energy Difference (meV/atom)")


# Murnaghan Equation of State Information
#           parameters = (E0, K0, K0', V0)

#   Diamond structure calculations
volumes_sim_diamond = np.power(lattice_parameters_diamond, 3) / number_of_atoms[0]
initial_parameters_diamond = (total_energies_strain_diamond[mid(total_energies_strain_diamond)], 1e11, 3.5, volumes_sim_diamond[mid(volumes_sim_diamond)])
fit_parameters_diamond = sp.fmin(square_differences, initial_parameters_diamond, args=(volumes_sim_diamond, total_energies_strain_diamond, murnaghan), maxiter=100000)
volumes_diamond = np.linspace(volumes_sim_diamond[0], volumes_sim_diamond[-1], num=100)
fit_total_energies_strain_diamond = murnaghan(fit_parameters_diamond, volumes_diamond)
# print(fit_parameters_diamond)

#   Diamond structure plot
figure_index += 1
fig7 = plt.figure(figure_index+1)
plt.plot(volumes_diamond, fit_total_energies_strain_diamond*6.242e18 / number_of_atoms[0])
plt.scatter(volumes_sim_diamond, total_energies_strain_diamond*6.242e18 / number_of_atoms[0])
plt.title(plot_titles[figure_index])
plt.xlabel(plot_xlabels[figure_index])
plt.ylabel(r'Total Energy (eV/atom)')


#     Beta-Sn structure calculations
volumes_sim_BetaSn = (celldm_3_BetaSn / number_of_atoms[1])*np.power(lattice_parameters_BetaSn, 3)
volumes_BetaSn = np.linspace(volumes_sim_BetaSn[0], volumes_sim_BetaSn[-1], num=100)
#
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

#   Beta-Sin structure plot
figure_index += 1
fig8 = plt.figure(figure_index+1)
plt.plot(volumes_BetaSn, fit_total_energies_strain_all_BetaSn*6.242e18 / number_of_atoms[1])
plt.scatter(volumes_sim_BetaSn, total_energies_strain_all_BetaSn*6.242e18 / number_of_atoms[1])
plt.plot(volumes_BetaSn, fit_total_energies_strain_shape_BetaSn * 6.242e18 / number_of_atoms[1])
plt.scatter(volumes_sim_BetaSn, total_energies_strain_shape_BetaSn * 6.242e18 / number_of_atoms[1])
plt.title(plot_titles[7])
plt.ylabel(r'Total Energy (eV/atom)')
plt.xlabel(r'Volume (m$^3$/atom)')

# Two-phase plot
figure_index += 1
fig9 = plt.figure(figure_index+1)
plt.plot(volumes_diamond/cubic_meters_per_cubic_angstrom, fit_total_energies_strain_diamond*6.242e18, label='Diamond')
plt.scatter(volumes_sim_diamond/cubic_meters_per_cubic_angstrom, total_energies_strain_diamond*6.242e18)
plt.plot(volumes_BetaSn/cubic_meters_per_cubic_angstrom, fit_total_energies_strain_all_BetaSn*6.242e18, label='BetaSn (cell\_dofree = all)')
plt.scatter(volumes_sim_BetaSn/cubic_meters_per_cubic_angstrom, total_energies_strain_all_BetaSn*6.242e18)
plt.plot(volumes_BetaSn/cubic_meters_per_cubic_angstrom, fit_total_energies_strain_shape_BetaSn*6.242e18, label='BetaSn (cell\_dofree = z)')
plt.scatter(volumes_sim_BetaSn/cubic_meters_per_cubic_angstrom, total_energies_strain_shape_BetaSn*6.242e18)
plt.yscale(r'linear')
plt.title(plot_titles[8])
plt.ylabel(r'Total Energy (eV/atom)')
plt.xlabel(r'Volume (\r{A}$^3$/atom)')
plt.legend()


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


# Transition volume
def transition_volume(a, p):
    e0 = a[0]
    k0 = a[1]
    k0prime = a[2]
    v0 = a[3]
    v = v0 * np.power((1 + p * (k0prime / k0)), (-1.0 / k0prime))
    return v
#     v = a[3]*np.power(((a[2]/a[1])*p+1),(-1/a[2]))
#     return a[0] + a[1]*a[3]*((1/(a[2]*((a[2]-1))))*np.power((v/a[3]), (1-a[3])) + (v/(a[2]*a[3])) - (1/(a[2]-1)))

    # return ((p[1]/p[2])*(np.power((v/p[3]),(-p[2]))-1))


fit_di = np.array([0,88.4e9,4.2,20.45e-30])
fit_be = np.array([.291/1.6e-19,106.1e9,4.6,15.34e-30])
div = 100000
vol_mash = np.concatenate((volumes_sim_diamond,volumes_sim_BetaSn))
volume = np.linspace(min(vol_mash), max(vol_mash), num=div)
pressure = np.linspace(-50*11.7e9 ,50*11.7e9,div)
diamond = enthal_murn(fit_parameters_diamond, pressure)
beta = enthal_murn(fit_parameters_shape_BetaSn, pressure)
# plt.figure()
# plt.plot(pressure, diamond)
# plt.plot(pressure, beta)
# plt.show()
matches = np.zeros(div)
index = np.arange(0,div,1)

for x in index:
    add = np.array([diamond[x]/beta[x]])
    matches[x] = add

diamond = diamond[np.logical_not(np.isnan(matches))]
beta = beta[np.logical_not(np.isnan(matches))]
pressure = pressure[np.logical_not(np.isnan(matches))]
matches = matches[np.logical_not(np.isnan(matches))]


# Find nearest neighbor atoms
def find_nearest(array, value):
    array = np.asarray(array)
    idx = np.argmin((np.abs(array - value)))
    return [array[idx], idx]


near = find_nearest(matches, 1)
print(near)
# tpressured = diamond[near[1]]
# tpressureb = beta[near[1]]
# tvolume = volume[near[1]]
# print(r'Transition pressure diamond:', tpressured/pascals_per_gigapascal, r'GPa')
# print(r'Transition pressure beta:', tpressureb/pascals_per_gigapascal, r'GPa')
# print(r'transition volume:', tvolume/cubic_meters_per_cubic_angstrom, r'Å^3')

tpressure = pressure[near[1]]
print('Transition pressure diamond:', tpressure/pascals_per_gigapascal, 'GPa')
tvol_diamond = transition_volume(fit_parameters_diamond, tpressure)
tvol_beta = transition_volume(fit_parameters_shape_BetaSn, tpressure)
print(tvol_diamond*1e30)
print(tvol_beta*1e30)
print((fit_parameters_shape_BetaSn[0]-fit_parameters_diamond[0])/1.6e-19)
print(fit_parameters_diamond)

plt.plot(volume/cubic_meters_per_cubic_angstrom, (murnaghan(fit_parameters_shape_BetaSn, tvol_beta)-(volume-tvol_beta)*tpressure)*6.242e18)

plt.show()