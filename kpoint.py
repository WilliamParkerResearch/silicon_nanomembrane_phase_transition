import matplotlib.pyplot as plt
import matplotlib as mpl
from QEData import *
from Information import *

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

figure_index = 2
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


plt.show()