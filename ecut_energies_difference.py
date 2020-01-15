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


indexes_ecut_diamond_BetaSn = matchmaker(cutoff_energies_diamond, cutoff_energies_BetaSn)
total_energies_ecut_diamond_matched = total_energies_kpoint_diamond[indexes_ecut_diamond_BetaSn[0]]
total_energies_ecut_BetaSn_matched = total_energies_kpoint_BetaSn[indexes_ecut_diamond_BetaSn[1]]
cutoff_energies_diamond_matched = cutoff_energies_diamond[indexes_ecut_diamond_BetaSn[0]]
cutoff_energies_BetaSn_matched = cutoff_energies_BetaSn[indexes_ecut_diamond_BetaSn[1]]

total_energies_ecut_diff = np.subtract(total_energies_ecut_BetaSn_matched, total_energies_ecut_diamond_matched)

# Finds appropriate y limits for plot
min_energy = np.sort(total_energies_ecut_diff)[0]
max_energy = np.sort(total_energies_ecut_diff)[-1]
one_interval = 0.01*(max_energy - min_energy)

y_minimum = min_energy - 5 * one_interval
y_maximum = max_energy + 5 * one_interval

#Plot
fig3 = plt.figure()
plt.plot(cutoff_energies_diamond_matched,total_energies_ecut_diff,label='energy difference')
plt.scatter(cutoff_energies_diamond_matched, total_energies_ecut_diff)
plt.title(r'Cutoff Energy Convergence of ' + structure_names[0] + ' and ' + structure_names[1] + ' ' +
               chemical_formula + ' Total Energy Difference')
plt.xlabel(r"Cutoff Energy (eV)")
plt.ylabel(r"Total Energy Difference (meV/atom)")
plt.xlim(cutoff_energies_diamond_matched[0], cutoff_energies_diamond_matched[-1])
plt.ylim(y_minimum, y_maximum)

plt.show()