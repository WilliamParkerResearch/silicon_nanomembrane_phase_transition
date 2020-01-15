import matplotlib.pyplot as plt
import matplotlib as mpl
from QEData import *


# Use TeX fonts
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.sans-serif'] = "cmr10"


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


# Finds appropriate y limits for plot
min_energy = np.sort(total_energies_ecut_diamond)[0]
max_energy = np.sort(total_energies_ecut_diamond)[-1]
one_interval = 0.01*(max_energy - min_energy)

y_minimum = min_energy - 5 * one_interval
y_maximum = max_energy + 5 * one_interval

#Plot
fig1 = plt.figure()
plt.plot(cutoff_energies_diamond, total_energies_ecut_diamond, label='total energy')
plt.scatter(cutoff_energies_diamond, total_energies_ecut_diamond)
plt.title(r'Cutoff Energy Convergence of Total Energy for ' + structure_names[0] + ' ' + chemical_formula)
plt.xlabel(r"Cutoff Energy (eV)")
plt.ylabel(r"Total Energy (meV/atom)")
plt.xlim(cutoff_energies_diamond[0], cutoff_energies_diamond[-1])
plt.ylim(y_minimum, y_maximum)
plt.axhline(total_energies_ecut_diamond[-1], color='r', label='convergence energy')
plt.axhline(total_energies_ecut_diamond[-1]-convergence_threshold, color='r', linestyle='--', label='convergence threshold')
plt.axhline(total_energies_ecut_diamond[-1]+convergence_threshold, color='r', linestyle='--')
plt.legend()

plt.show()