import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib as mpl
from QEData import *

# Get plot parameters
from format_charts import *

# Decide output form
save_plot = True
plot_filename = 'Si.I41amd.PBE.bands.pdf'
show_plot = False


num_disc = len(disc_points_betasn)
index_array = np.array([])
last_index = int(len(kpoint_indexes_betasn)-1)

for i in disc_points_betasn:
    index = int(np.where(kpoint_indexes_betasn == i)[0][0])
    index_array = np.append(index_array, index)

index_array = np.append(index_array, last_index)


for i in range(n_bands_betasn):
    for j in range(num_disc):
        if j == 0:
            exec(f'eigenarray_{i+1}_{j} = bands_eigenvalues_betasn_{i+1}[int(index_array[{j}]):int(index_array[{j+1}]+1)]')
            exec(f'karray_{j} = kpoint_indexes_betasn[int(index_array[{j}]):int(index_array[{j+1}]+1)]')

        else:
            exec(f'eigenarray_{i+1}_{j} = bands_eigenvalues_betasn_{i+1}[int(index_array[{j}]+1):int(index_array[{j+1}]+1)]')
            exec(f'karray_{j} = kpoint_indexes_betasn[int(index_array[{j}]+1):int(index_array[{j+1}]+1)]')


disc_index_array = []
disc_last_index = int(len(ticks_betasn)-1)

for i in disc_points_betasn:
    disc_index = int(np.where(ticks_betasn == i)[0][0])
    disc_index_array = np.append(disc_index_array, disc_index)

disc_index_array = np.append(disc_index_array, disc_last_index)

ratio_array = np.diff(disc_index_array)

for i in ratio_array:
    if (i != ratio_array[0]):
        location = int(np.where(ratio_array == i)[0][0])
        ratio_array[location] = ratio_array[location] - 1               #to account for the fact that the discontinuty point is couted twice

ratio_array = np.append(ratio_array,1)

min_eigenvalue = np.sort(bands_eigenvalues_betasn_1)[0] - fermi_energy_betasn
exec(f'max_eigenvalue = np.sort(bands_eigenvalues_betasn_{n_bands_betasn})[-1] - {fermi_energy_betasn}')


fig = plt.figure()
# fig.suptitle(r'$\vec{k}$-point and Density of State Electrical Energy Bands for ' + structure_names[1] + ' ' + chemical_formula, )
gs = GridSpec(1, len(disc_points_betasn)+1, width_ratios=ratio_array)

for i in range(num_disc):
    if (i == 0):
        start = np.where(ticks_betasn == disc_points_betasn[i])[0][0]
        end = np.where(ticks_betasn == disc_points_betasn[i+1])[0][0]+1
    elif (i == num_disc - 1):
        start = np.where(ticks_betasn == disc_points_betasn[i])[0][0] +1
        end = len(ticks_betasn)
    else:
        start = np.where(ticks_betasn == disc_points_betasn[i])[0][0] +1
        end = np.where(ticks_betasn == disc_points_betasn[i + 1])[0][0] +1
    for j in range(n_bands_betasn):
        if (i == 0):
            exec(f'ax{i} = fig.add_subplot(gs[{i}])')
            exec(f'ax{i}.plot(karray_{i}, eigenarray_{j+1}_{i} - fermi_energy_betasn, color=band_color)')
            exec(f'ax{i}.set_xticks(ticks_betasn[start:end])')
            exec(f'ax{i}.set_xticklabels(ticklabels_betasn[start:end])')
            exec(f'ax{i}.set_xlim([disc_points_betasn[{i}],disc_points_betasn[{i+1}]])')
            exec(f'ax{i}.set_ylim([min_eigenvalue,max_eigenvalue])')
            #exec(f'ax{i}.set_ylabel(r\'Energy (eV)\')')
            ax0.set_ylabel(r'$\varepsilon_i(\vec{k})$')
            ax0.set_xlabel(r'$\vec{k}$')
        elif (i == num_disc -1):
            exec(f'ax{i} = fig.add_subplot(gs[{i}], sharey = ax0)')
            exec(f'ax{i}.plot(karray_{i}, eigenarray_{j+1}_{i} - fermi_energy_betasn, color=band_color)')
            exec(f'ax{i}.tick_params(labelleft=False)')
            exec(f'ax{i}.set_xticks(ticks_betasn[start:end])')
            exec(f'ax{i}.set_xticklabels(ticklabels_betasn[start:end])')
            exec(f'ax{i}.set_xlim([disc_points_betasn[{i}],kpoint_indexes_betasn[-1]])')
        else:
            exec(f'ax{i} = fig.add_subplot(gs[{i}], sharey = ax0)')
            exec(f'ax{i}.plot(karray_{i}, eigenarray_{j+1}_{i} - fermi_energy_betasn, color=band_color)')
            exec(f'ax{i}.tick_params(labelleft=False)')
            exec(f'ax{i}.set_xticks(ticks_betasn[start:end])')
            exec(f'ax{i}.set_xticklabels(ticklabels_betasn[start:end])')
            exec(f'ax{i}.set_xlim([disc_points_betasn[{i}],disc_points_betasn[{i+1}]])')
        exec(f'ax{i}.axhline(y=0, color=fermi_level_linecolor, linestyle=fermi_level_linestyle, linewidth=fermi_level_linewidth)')
ax_last = fig.add_subplot(gs[num_disc], sharey = ax0)
ax_last.tick_params(labelleft=False)
ax_last.plot(density_betasn,dos_energies_betasn - fermi_energy_betasn, color=dos_curve_color)
ax_last.fill(density_betasn, dos_energies_betasn - fermi_energy_betasn, color=dos_fill_color, alpha=dos_opacity)
# ax_last.set_xlabel(r'Electric Density (E)')
ax_last.set_xlabel(r'$g(\varepsilon)$')
plt.subplots_adjust(wspace=0.05)

if save_plot:
    plt.savefig(plot_filename)

if show_plot:
    plt.show()
