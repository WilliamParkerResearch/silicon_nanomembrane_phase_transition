import numpy as np
from ReferenceFiles.vasp_python_converter import vasp_data_modifier, bulk_cell_vstacker
from ReferenceFiles.pair_distribution import *
from ReferenceFiles.plot_formating import *
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

exchange_correlation = 'PBE'
phase = 'diamond'
# phase = 'betasn'

if phase == 'diamond':
    prefix = 'Si.Fd-3m_'
    density = 2300
    rgbcode = rgbcode_diamond
if phase == 'betasn':
    prefix = 'Si.I4_1amd_'
    density = 3370
    rgbcode = rgbcode_betasn


verbose = False
save_plot = False
minimum_distance = 1.5
maximum_distance = 6
# minimum_distance = 2.35
# maximum_distance = 2.6
l_n_ml = 0
h_n_ml = 8
# np.arange(l_n_ml, h_n_ml+1)

fig = plt.figure(figsize=fs_m_12)
layers = np.array([0,1,2,h_n_ml]) #np.arange(l_n_ml, h_n_ml+1) #
n_entries = len(layers)

k = 0


for i in layers:
    output_file = 'DataFolder'+'/'+exchange_correlation+'/'+'output_files'+'/'+prefix+str(i)+'ML.scf.out'
    # lattice_parameter = get_celldm(output_file, verbose=verbose)
    # lattice_parameter *= constants.value('Bohr radius')/constants.angstrom
    # atomic_positions = lattice_parameter * multislab_builder(get_positions(output_file, verbose=verbose), 3, 3)
    # atomic_positions = lattice_parameter * get_positions(output_file)
    lattice_parameter = vasp_data_modifier(i,phase)[1]
    atomic_positions = vasp_data_modifier(i,phase)[0]
    lattice_parameter0 = vasp_data_modifier(0,phase)[1]
    atomic_positions0 = vasp_data_modifier(0,phase)[0]
    if i == 0:
        atomic_positions = bulk_cell_vstacker(h_n_ml,phase)


    if maximum_distance <= 0.05:
        maximum_distance = np.ceil(lattice_parameter)

    if i == 0:

        rgbacode = adjust_lightness(rgbcode,0.5)
        plot_pdf(atomic_positions, density, zorder=n_entries+1, minimum_distance=minimum_distance,maximum_distance=maximum_distance,color=rgbacode,linewidth=universal_linewidth+0.3,dashes=(5,5))
    else:
        minalpha = 2.3
        topalpha = 0.7
        pdiv = (1/(n_entries))
        alpha = minalpha+k*pdiv*(topalpha-minalpha)
        rgbacode = adjust_lightness(rgbcode,alpha)
        plot_pdf(atomic_positions, density, zorder=k+1,minimum_distance=minimum_distance,maximum_distance=maximum_distance,linestyle='solid',color=rgbacode,linewidth=universal_linewidth)
    k+=1


plt.xlabel(r'$r ({\rm \AA})$',fontsize=axis_fontsize)
plt.ylabel(r'$g(r) $',fontsize=axis_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.xlim([minimum_distance, maximum_distance])

lines = [Line2D([0], [0], color=adjust_lightness(rgbcode,1.4),linewidth=universal_linewidth),
         Line2D([0], [0], color=adjust_lightness(rgbcode,0.6),linestyle='dashed',linewidth=universal_linewidth)]
labels = [r'$N_{ML}$','Bulk']

plt.legend(lines,labels,loc='upper right')
plot(False,'histogram_'+phase)