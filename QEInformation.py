import numpy as np
import scipy.optimize as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
from QEData import *

# Use TeX fonts
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.sans-serif'] = "cmr10"

#constants
convergence_threshold = 7.35e-5     #this is = 1meV

#definitions
def murnaghan(p, v):
    kk = p[2]-1.0
    return (p[0] + (p[1]*p[3]*(((1.0/(p[2]*kk))*np.power((v/p[3]), (-kk)))+(v/(p[2]*p[3]))-(1.0/kk))))


def square_differences(p, x, y, f):
    return np.power(f(p, x) - y, 2).sum()


def mid(x):
    return int((len(x)-1)/2)


def matchmaker(a,b):

    list_a = list(a)
    list_b = list(b)

    matches = np.intersect1d(a,b)
    index_a = []
    for x in matches:
        index_a.append(list_a.index(x))
    index_b = []
    for x in matches:
        index_b.append(list_b.index(x))
    return [index_a, index_b]


#ecut information
    #diamond
fig1 = plt.figure(1)
plt.plot(cutoff_energies_diamond, total_energies_ecut_diamond, label='Total energy')
plt.scatter(cutoff_energies_diamond, total_energies_ecut_diamond)
plt.title(r"Cutoff Energy vs Total Energy for Diamond Structure of Si")
plt.xlabel(r"Cutoff Energy (meV)")
plt.ylabel(r"Total Energy (meV/atom)")
plt.axhline(total_energies_ecut_diamond[-1], color='r', label='convergence energy')
plt.axhline(total_energies_ecut_diamond[-1]-convergence_threshold, color='r', linestyle='--', label='convergence threshold')
plt.axhline(total_energies_ecut_diamond[-1]+convergence_threshold, color='r', linestyle='--')
plt.legend()
# plt.axis([0, ecut[-1],etot[-1]-10*convergence_threshold, etot[-1]+10*convergence_threshold])
    #BetaSn
fig2 = plt.figure(2)
plt.plot(cutoff_energies_BetaSn, total_energies_ecut_BetaSn, label='Total energy')
plt.scatter(cutoff_energies_BetaSn, total_energies_ecut_BetaSn)
plt.title("Cutoff Energy vs Total Energy for BetaSn Structure of Si")
plt.xlabel("Cutoff Energy (meV)")
plt.ylabel("Total Energy (meV/atom)")
plt.axhline(total_energies_ecut_BetaSn[-1], color='r', label='convergence energy')
plt.axhline(total_energies_ecut_BetaSn[-1]-convergence_threshold, color='r', linestyle='--', label='convergence threshold')
plt.axhline(total_energies_ecut_BetaSn[-1]+convergence_threshold, color='r', linestyle='--')
plt.legend()
# plt.axis([0, ecut[-1],etot[-1]-10*convergence_threshold, etot[-1]+10*convergence_threshold])

    #Difference Between Diamond and BetaSn
indexes_ecut_diamond_BetaSn = matchmaker(cutoff_energies_diamond, cutoff_energies_BetaSn)
total_energies_ecut_diamond_matched = total_energies_kpoint_diamond[indexes_ecut_diamond_BetaSn[0]]
total_energies_ecut_BetaSn_matched = total_energies_kpoint_BetaSn[indexes_ecut_diamond_BetaSn[1]]
cutoff_energies_diamond_matched = cutoff_energies_diamond[indexes_ecut_diamond_BetaSn[0]]
cutoff_energies_BetaSn_matched = cutoff_energies_BetaSn[indexes_ecut_diamond_BetaSn[1]]

total_energies_ecut_diff = np.subtract(total_energies_ecut_BetaSn_matched, total_energies_ecut_diamond_matched)

fig3 = plt.figure(3)
plt.plot(cutoff_energies_diamond_matched,total_energies_ecut_diff,label='energy difference')
plt.scatter(cutoff_energies_diamond_matched, total_energies_ecut_diff)
plt.title(r"Cutoff Energy vs Total Energy Differences Between Diamond Form and Beta-tin Form")
plt.xlabel(r"Cutoff Energy (meV)")
plt.ylabel(r"Total Energy Difference (meV/atom)")


#Kpoint Information
    #diamond
fig4 = plt.figure(4)
plt.plot(kpoints_diamond,total_energies_kpoint_diamond, label='total energy')
plt.scatter(kpoints_diamond, total_energies_kpoint_diamond)
plt.title(r"Kpoints vs Total Energy of Diamond Structure of Si")
plt.xlabel(r"Number of Kpoints")
plt.ylabel(r"Total Energy (meV/atom)")
plt.axhline(total_energies_kpoint_diamond[-1], color='r', label='convergence energy')
plt.axhline(total_energies_kpoint_diamond[-1]-convergence_threshold, color='r', linestyle='--', label='convergence threshold')
plt.axhline(total_energies_kpoint_diamond[-1]+convergence_threshold, color='r', linestyle='--')
plt.axis([0, kpoints_diamond[-1],total_energies_kpoint_diamond[-1]-10*convergence_threshold, total_energies_kpoint_diamond[-1]+10*convergence_threshold])
plt.legend()

    #BetaSn
fig5 = plt.figure(5)
plt.plot(kpoints_BetaSn,total_energies_kpoint_BetaSn, label='total energy')
plt.scatter(kpoints_BetaSn, total_energies_kpoint_BetaSn)
plt.title(r"Kpoints vs Total Energy of BetaSn Structure of Si")
plt.xlabel(r"Number of Kpoints")
plt.ylabel(r"Total Energy (meV/atom)")
plt.axhline(total_energies_kpoint_BetaSn[-1], color='r', label='convergence energy')
plt.axhline(total_energies_kpoint_BetaSn[-1]-convergence_threshold, color='r', linestyle='--', label='convergence threshold')
plt.axhline(total_energies_kpoint_BetaSn[-1]+convergence_threshold, color='r', linestyle='--')
plt.axis([0, kpoints_BetaSn[-1],total_energies_kpoint_BetaSn[-1]-10*convergence_threshold, total_energies_kpoint_BetaSn[-1]+10*convergence_threshold])
plt.legend()

    #Difference Between Diamond and BetaSn
indexes_kpoint_diamond_BetaSn = matchmaker(kpoints_diamond, kpoints_BetaSn)
total_energies_kpoint_diamond_matched = total_energies_kpoint_diamond[indexes_kpoint_diamond_BetaSn[0]]
total_energies_kpoint_BetaSn_matched = total_energies_kpoint_BetaSn[indexes_kpoint_diamond_BetaSn[1]]
kpoints_diamond_matched = kpoints_diamond[indexes_kpoint_diamond_BetaSn[0]]
kpoints_BetaSn_matched = kpoints_BetaSn[indexes_kpoint_diamond_BetaSn[1]]

total_energies_kpoint_diff = np.subtract(total_energies_kpoint_BetaSn_matched, total_energies_kpoint_diamond_matched)

fig6 = plt.figure(6)
plt.plot(kpoints_diamond_matched,total_energies_kpoint_diff)
plt.scatter(kpoints_diamond_matched, total_energies_kpoint_diff)
plt.title(r"Kpoints vs Total Energy Differences Between Diamond Form and Beta-tin Form")
plt.xlabel(r"Number of Kpoints")
plt.ylabel(r"Total Energy Difference (meV/atom)")


#Murnaghan Equation of State Information
# parameters = (E0, K0, K0', V0)

    #diamond
volumes_sim_diamond = (1/8)*np.power(lattice_parameters_diamond, 3)
initial_parameters_diamond = (total_energies_strain_diamond[mid(total_energies_strain_diamond)], 1e11, 3.5, volumes_sim_diamond[mid(volumes_sim_diamond)])
fit_parameters_diamond = sp.fmin(square_differences, initial_parameters_diamond, args=(volumes_sim_diamond, total_energies_strain_diamond, murnaghan), maxiter=100000)
volumes_diamond = np.linspace(volumes_sim_diamond[0], volumes_sim_diamond[-1], num=100)
fit_total_energies_strain_diamond = murnaghan(fit_parameters_diamond, volumes_diamond)

# print(fit_parameters_diamond)

fig7 = plt.figure(7)
plt.plot(volumes_diamond, fit_total_energies_strain_diamond*6.242e18/8)
plt.scatter(volumes_sim_diamond, total_energies_strain_diamond*6.242e18/8)
plt.title(r'Equation of State of Diamond Structure')
plt.ylabel(r'Total Energy (eV/atom)')
plt.xlabel(r'Volume (m$^3$/atom)')

#     #BetaSn
volumes_sim_BetaSn = (celldm_3_BetaSn/4)*np.power(lattice_parameters_BetaSn, 3)
volumes_BetaSn = np.linspace(volumes_sim_BetaSn[0], volumes_sim_BetaSn[-1], num=100)
#
          #cell_dofree ='all'
initial_parameters_all_BetaSn = (total_energies_strain_all_BetaSn[mid(total_energies_strain_all_BetaSn)], 1e11, 3.5, volumes_sim_BetaSn[mid(volumes_sim_BetaSn)])
fit_parameters_all_BetaSn = sp.fmin(square_differences, initial_parameters_all_BetaSn, args=(volumes_sim_BetaSn, total_energies_strain_all_BetaSn, murnaghan), maxiter=100000)
fit_total_energies_strain_all_BetaSn = murnaghan(fit_parameters_all_BetaSn, volumes_BetaSn)

# print(fit_parameters_all_BetaSn)

        #cell_dofree ='z'
initial_parameters_z_BetaSn = (total_energies_strain_z_BetaSn[mid(total_energies_strain_z_BetaSn)], 1e11, 3.5, volumes_sim_BetaSn[mid(volumes_sim_BetaSn)])
fit_parameters_z_BetaSn = sp.fmin(square_differences, initial_parameters_z_BetaSn, args=(volumes_sim_BetaSn, total_energies_strain_z_BetaSn, murnaghan), maxiter=100000)
fit_total_energies_strain_z_BetaSn = murnaghan(fit_parameters_z_BetaSn, volumes_BetaSn)

# print(fit_parameters_z_BetaSn)

fig8 = plt.figure(8)
plt.plot(volumes_BetaSn, fit_total_energies_strain_all_BetaSn*6.242e18/4)
plt.scatter(volumes_sim_BetaSn, total_energies_strain_all_BetaSn*6.242e18/4)
plt.plot(volumes_BetaSn, fit_total_energies_strain_z_BetaSn * 6.242e18/4)
plt.scatter(volumes_sim_BetaSn, total_energies_strain_z_BetaSn * 6.242e18/4)
plt.title(r'Equation of State of BetaSn Structure')
plt.ylabel(r'Total Energy (eV/atom)')
plt.xlabel(r'Volume (m$^3$/atom)')


            #multiple plot
fig9 = plt.figure(9)
plt.plot(volumes_diamond/1e-30, fit_total_energies_strain_diamond*6.242e18, label='Diamond')
plt.scatter(volumes_sim_diamond/1e-30, total_energies_strain_diamond*6.242e18)
plt.plot(volumes_BetaSn/1e-30, fit_total_energies_strain_all_BetaSn*6.242e18, label='BetaSn (cell\_dofree = all)')
plt.scatter(volumes_sim_BetaSn/1e-30, total_energies_strain_all_BetaSn*6.242e18)
plt.plot(volumes_BetaSn/1e-30, fit_total_energies_strain_z_BetaSn*6.242e18, label='BetaSn (cell\_dofree = z)')
plt.scatter(volumes_sim_BetaSn/1e-30, total_energies_strain_z_BetaSn*6.242e18)
plt.yscale(r'linear')
plt.title(r'Equation of State')
plt.ylabel(r'Total Energy (eV/atom)')
plt.xlabel(r'Volume (\r{A}$^3$/atom)')
plt.legend()

#transition properties

# parameters = (E0, K0, K0', V0)
def murnaghan_prime(p, v):
    return (p[1]/p[2])*(1-np.power((v/p[3]),-p[2]))
div = 1000000
vol_mash = np.concatenate((volumes_sim_diamond,volumes_sim_BetaSn))
volume  = np.linspace(0.25*min(vol_mash), 1.75*max(vol_mash), num=div)
diamond = murnaghan_prime(fit_parameters_diamond, volume)
beta = murnaghan_prime(fit_parameters_z_BetaSn, volume)

matches = np.zeros(div)
index = np.arange(0,div,1)

for x in index:
    add = np.array([diamond[x]/beta[x]])
    matches[x] = add

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return [array[idx], idx]

near = find_nearest(matches, 1)
print(near)
tpressured = diamond[near[1]]
tpressureb = beta[near[1]]
tvolume = volume[near[1]]
print(r'Transition pressure diamond:', tpressured/1e9, r'GPa')
print(r'Transition pressure beta:', tpressureb/1e9, r'GPa')
print(r'transition volume:', tvolume/1e-30, r'Å^3')

plt.show()