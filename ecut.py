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


plt.show()