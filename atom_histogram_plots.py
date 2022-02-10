import numpy as np
from ReferenceFiles.plot_formating import *
from ReferenceFiles.histogram_functions import plot_histogram, plot_histogram_differences
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


show_plot = True
exchange_correlation = 'PBE'
phase = 'diamond'
# phase = 'betasn'
proj = 'n'
l_n_ml = 0
h_n_ml = 8
# layers = np.array([0,5]) #np.arange(l_n_ml, h_n_ml+1)
layers = np.array([0,1,3,8]) #np.arange(l_n_ml, h_n_ml+1)
ind_plot=True
ind_atoms = ['cut',0]
max_len = 6
min_len = 0.3

if phase == 'diamond':
    rgbcode = rgbcode_diamond
if phase == 'betasn':
    rgbcode = rgbcode_betasn





if len(layers)> 5:
    lines = [Line2D([0], [0], color=adjust_lightness(rgbcode, 1.4), linewidth=universal_linewidth),
             Line2D([0], [0], color=adjust_lightness(rgbcode, 0.6), linestyle='dashed', linewidth=universal_linewidth)]
    labels = [r'$N_{ML}$', 'Bulk']



# # ########################### Regular histogram plot ################################
# fig = plt.figure(figsize=fs_m_12)
# plot_histogram(layers,phase)
#
# if len(layers) > 5:
#     plt.legend(lines,labels,loc='upper right')
# else:
#     plt.legend()
#
# plot(show_plot,'histogram_'+phase)
#
#
# # ######################## Histogram differences plot #######################################
# fig = plt.figure(figsize=fs_m_12)
# plot_histogram_differences(layers,phase)
#
# if len(layers) > 5:
#     plt.legend(lines,labels,loc='upper right')
# else:
#     plt.legend()
#
# plot(False,'histogram_differences_'+phase)


########################### Histogram joint plot ####################################

fig = plt.figure(figsize=fs_m_23)
# plt.ylabel(r'${g(r)}_{{N}_{ML}} $',fontsize=axis_fontsize)
plot_histogram_differences(layers,phase,proj=proj,maximum_distance=max_len,minimum_distance=min_len,individual_atom_plot=ind_plot,plot_n_atoms=ind_atoms)








if len(layers) > 5:
    plt.legend(lines,labels,loc='upper right')
else:
    plt.legend()

plot(True,'histogram_subplot_'+phase)
