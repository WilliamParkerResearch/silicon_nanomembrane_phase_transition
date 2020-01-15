import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib as mpl
from QEData import *



# Use TeX fonts
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.sans-serif'] = "cmr10"

num_disc = len(disc_points_diamond)
index_array = np.array([])
last_index = int(len(kpoint_indexes_diamond)-1)

for i in disc_points_diamond:
    index = int(np.where(kpoint_indexes_diamond == i)[0][0])
    index_array = np.append(index_array, index)

index_array = np.append(index_array, last_index)


for i in range(n_bands_diamond):
    for j in range(num_disc):
        if j == 0:
            exec(f'eigenarray_{i+1}_{j} = bands_eigenvalues_diamond_{i+1}[int(index_array[{j}]):int(index_array[{j+1}]+1)]')
            exec(f'karray_{j} = kpoint_indexes_diamond[int(index_array[{j}]):int(index_array[{j+1}]+1)]')

        else:
            exec(f'eigenarray_{i+1}_{j} = bands_eigenvalues_diamond_{i+1}[int(index_array[{j}]+1):int(index_array[{j+1}]+1)]')
            exec(f'karray_{j} = kpoint_indexes_diamond[int(index_array[{j}]+1):int(index_array[{j+1}]+1)]')



disc_index_array = []
disc_last_index = int(len(ticks_diamond)-1)

for i in disc_points_diamond:
    disc_index = int(np.where(ticks_diamond == i)[0][0])
    disc_index_array = np.append(disc_index_array, disc_index)

disc_index_array = np.append(disc_index_array, disc_last_index)

ratio_array = np.diff(disc_index_array)

for i in ratio_array:
    if (i != ratio_array[0]):
        location = int(np.where(ratio_array == i)[0][0])
        ratio_array[location] = ratio_array[location] - 1               #to account for the fact that the discontinuty point is couted twice

ratio_array = np.append(ratio_array,1)

min_eigenvalue = np.sort(bands_eigenvalues_diamond_1)[0] - fermi_energy_diamond
exec(f'max_eigenvalue = np.sort(bands_eigenvalues_diamond_{n_bands_diamond})[-1] - {fermi_energy_diamond}')


fig = plt.figure()
plt.suptitle(r'$\vec{k}$-point and Density of State Electrical Energy Bands for ' + structure_names[0] + ' ' + chemical_formula)
gs = GridSpec(1, len(disc_points_diamond)+1, width_ratios=ratio_array)

for i in range(num_disc):
    if (i == 0):
        start = np.where(ticks_diamond == disc_points_diamond[i])[0][0]
        end = np.where(ticks_diamond == disc_points_diamond[i+1])[0][0]+1
    elif (i == num_disc - 1):
        start = np.where(ticks_diamond == disc_points_diamond[i])[0][0] +1
        end = len(ticks_diamond)
    else:
        start = np.where(ticks_diamond == disc_points_diamond[i])[0][0] +1
        end = np.where(ticks_diamond == disc_points_diamond[i + 1])[0][0]

    for j in range(n_bands_diamond):
        if (i == 0):
            exec(f'ax{i} = fig.add_subplot(gs[{i}])')
            exec(f'ax{i}.plot(karray_{i}, eigenarray_{j+1}_{i} - fermi_energy_diamond)')
            exec(f'ax{i}.set_xticks(ticks_diamond[start:end])')
            exec(f'ax{i}.set_xticklabels(ticklabels_diamond[start:end])')
            exec(f'ax{i}.set_xlim([disc_points_diamond[{i}],disc_points_diamond[{i+1}]])')
            exec(f'ax{i}.set_ylim([min_eigenvalue,max_eigenvalue])')
            exec(f'ax{i}.set_ylabel(r\'Energy (eV)\')')
            ax0.set_xlabel(r'$\vec{k}$ points')
        elif (i == num_disc -1):
            exec(f'ax{i} = fig.add_subplot(gs[{i}], sharey = ax0)')
            exec(f'ax{i}.plot(karray_{i}, eigenarray_{j+1}_{i} - fermi_energy_diamond)')
            exec(f'ax{i}.tick_params(labelleft=False)')
            exec(f'ax{i}.set_xticks(ticks_diamond[start:end])')
            exec(f'ax{i}.set_xticklabels(ticklabels_diamond[start:end])')
            exec(f'ax{i}.set_xlim([disc_points_diamond[{i}],kpoint_indexes_diamond[-1]])')
        else:
            exec(f'ax{i} = fig.add_subplot(gs[{i}], sharey = ax0)')
            exec(f'ax{i}.plot(karray_{i}, eigenarray_{j+1}_{i} - fermi_energy_diamond)')
            exec(f'ax{i}.tick_params(labelleft=False)')
            exec(f'ax{i}.set_xticks(ticks_diamond[start:end])')
            exec(f'ax{i}.set_xticklabels(ticklabels_diamond[start:end])')
            exec(f'ax{i}.set_xlim([disc_points_diamond[{i}],disc_points_diamond[{i+1}]])')
        exec(f'ax{i}.axhline(y=0, color=\'r\', linestyle=\'--\', linewidth=0.5)')
ax_last = fig.add_subplot(gs[num_disc], sharey = ax0)
ax_last.tick_params(labelleft=False)
ax_last.plot(density_diamond,dos_energies_diamond - fermi_energy_diamond)
ax_last.set_xlabel(r'Electric Density (E)')
plt.subplots_adjust(wspace=0.05)
plt.show()

print(ratio_array)