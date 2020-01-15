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
min_energy = np.sort(total_energies_kpoint_BetaSn)[0]
max_energy = np.sort(total_energies_kpoint_BetaSn)[-1]
one_interval = 0.01*(max_energy - min_energy)

y_minimum = min_energy - 5 * one_interval
y_maximum = max_energy + 5 * one_interval

#Plot
fig5 = plt.figure()
plt.plot(kpoints_BetaSn,total_energies_kpoint_BetaSn, label='total energy')
plt.scatter(kpoints_BetaSn, total_energies_kpoint_BetaSn)
plt.title(r'$\vec{k}$-point Mesh Convergence of Total Energy for ' + structure_names[1] + ' ' + chemical_formula)
plt.xlabel(r"Number of $\vec{k}$ points in each direction")
plt.ylabel(r"Total Energy (meV/atom)")
plt.axhline(total_energies_kpoint_BetaSn[-1], color='r', label='convergence energy')
plt.axhline(total_energies_kpoint_BetaSn[-1]-convergence_threshold, color='r', linestyle='--', label='convergence threshold')
plt.axhline(total_energies_kpoint_BetaSn[-1]+convergence_threshold, color='r', linestyle='--')
plt.xlim(kpoints_BetaSn[0], kpoints_BetaSn[-1])
plt.ylim(y_minimum, y_maximum)
plt.legend()


plt.show()