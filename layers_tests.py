import numpy as np
import matplotlib.pyplot as plt
from eos_information import eos_properties
from ReferenceFiles.cohesive_energy_values import *
from ReferenceFiles.unit_conversions import *
from scipy.optimize import curve_fit, fmin
from ReferenceFiles.FunctionDefinitions import mid, square_differences
from ReferenceFiles.plot_formating import *
from matplotlib.lines import Line2D
from ReferenceFiles.equations_of_state import pressure_from_energy_equation_of_state
from ReferenceFiles.vasp_python_converter import vasp_data_modifier, bulk_cell_vstacker, ind_atom_distances, ztest
from ReferenceFiles.pair_distribution import format_pdf, plot_pdf
from scipy.signal import find_peaks


XC = 'PBE'
fold = 'eos'
showplot=True
der_acc = 1e-40
# showplot=False
#
#Use TeX fonts
# mpl.rcParams['text.usetex'] = True
# mpl.rcParams['font.sans-serif'] = "cmr10"
# mpl.rcParams['font.family'] = 'serif'
# mpl.rcParams['axes.formatter.useoffset'] = False

adj_l = 1.2
linewidth = universal_linewidth

n_layers = np.arange(0,9,1)
# n_layers = np.array([0,1,2])
# n_layers = np.array([0,1,2,7])

transition_pressures = np.array([])
transition_volumes_diamond = np.array([])
transition_volumes_betasn = np.array([])
volumes0_diamond = np.array([])
volumes0_betasn = np.array([])
energy0_diamond = np.array([])
energy0_betasn = np.array([])
bulk_modulus_diamond = np.array([])
bulk_modulus_betasn = np.array([])
bulk_modulus_prime_diamond = np.array([])
bulk_modulus_prime_betasn = np.array([])
cell_a_diamond = np.array([])
cell_c_diamond = np.array([])
cell_a_betasn = np.array([])
cell_c_betasn = np.array([])




for i in n_layers:

    if i == 0:
        t_pres0, vol_d0, vol_b0, fitp_d0, fitp_b0,ad0,ab0,cd0,cb0,_,_ = eos_properties(i,exchange_correlation=XC,eos_folder=fold)
        vol0_d0 = fitp_d0[3]
        energy0_d0 = fitp_d0[0]
        bulk_mod0_d0 = fitp_d0[1]
        bulk_modp0_d0 = fitp_d0[2]
        vol0_b0 = fitp_b0[3]
        energy0_b0 = fitp_b0[0]
        bulk_mod0_b0 = fitp_b0[1]
        bulk_modp0_b0 = fitp_b0[2]
        max_bond_d0 = max_bond_d
        max_bond_b0 = max_bond_b
        min_bond_b0 = min_bond_b

        dP_dV_d0 = np.diff(pressure_from_energy_equation_of_state(fitp_d0,[fitp_d0[3],fitp_d0[3]+der_acc]))/der_acc
        dP_dV_b0 = np.diff(pressure_from_energy_equation_of_state(fitp_b0,[fitp_b0[3],fitp_b0[3]+der_acc]))/der_acc

    else:
        t_pres, vol_d, vol_b, fitp_d, fitp_b,ad,ab,cd,cb,_,_ = eos_properties(i,exchange_correlation=XC,eos_folder=fold)
        transition_pressures = np.append(transition_pressures,t_pres)
        transition_volumes_diamond = np.append(transition_volumes_diamond,vol_d)
        transition_volumes_betasn = np.append(transition_volumes_betasn,vol_b)
        volumes0_diamond = np.append(volumes0_diamond,fitp_d[3])
        volumes0_betasn = np.append(volumes0_betasn,fitp_b[3])
        energy0_diamond = np.append(energy0_diamond,fitp_d[0])
        energy0_betasn = np.append(energy0_betasn,fitp_b[0])
        bulk_modulus_diamond = np.append(bulk_modulus_diamond,fitp_d[1])
        bulk_modulus_betasn = np.append(bulk_modulus_betasn,fitp_b[1])
        bulk_modulus_prime_diamond = np.append(bulk_modulus_prime_diamond,fitp_d[2])
        bulk_modulus_prime_betasn = np.append(bulk_modulus_prime_betasn,fitp_b[2])
        cell_a_diamond = np.append(cell_a_diamond,ad)
        cell_a_betasn = np.append(cell_a_betasn,ab)
        cell_c_diamond = np.append(cell_c_diamond, cd)
        cell_c_betasn = np.append(cell_c_betasn, cb)
        bonds_diamond = np.append(bonds_diamond, max_bond_d)
        bonds_betasn_max = np.append(bonds_betasn_max, max_bond_b)
        bonds_betasn_min = np.append(bonds_betasn_min, min_bond_b)


cohesive_energy0_diamond = -((energy0_diamond)-atom_energy_pbe)
cohesive_energy0_d0 = -((energy0_d0)-atom_energy_pbe)
cohesive_energy0_betasn = -((energy0_betasn)-atom_energy_pbe)
cohesive_energy0_b0 = -((energy0_b0)-atom_energy_pbe)





#######Diamond Plots
# ####Cohesive Energies
# fig = plt.figure(figsize=fs_s_11)
# plt.plot(n_layers[1:],cohesive_energy0_diamond*eV_per_joule,zorder=2,color=rgbcode_diamond_m,linewidth=linewidth)
# plt.scatter(n_layers[1:],cohesive_energy0_diamond*eV_per_joule,zorder=3,color=rgbcode_diamond_m,linewidth=linewidth, s=marker_size,marker=mark_d)
# plt.hlines(y=cohesive_energy0_d0*eV_per_joule,xmin=n_layers[1],xmax=n_layers[-1],zorder=1,color=adjust_lightness(rgbcode_diamond,adj_l),linewidth=linewidth,linestyle='dashed')
# plt.xlabel('N-ML',fontsize=axis_fontsize)
# plt.ylabel(r'$\mathit{E}$ (eV/atom)',fontsize=axis_fontsize)
# plt.xticks(fontsize=tick_fontsize)
# plt.yticks(fontsize=tick_fontsize)
#
# lines = [Line2D([0], [0], color=rgbcode_diamond_m,marker=mark_d,linewidth=universal_linewidth),
#          Line2D([0], [0], color=adjust_lightness(rgbcode_diamond,adj_l),linestyle='dashed',linewidth=universal_linewidth)]
# labels = ['N-ML','Bulk']
#
# plt.legend(lines,labels,loc='upper right')
# plot(showplot,'nml_cenergies_diamond')
#
# ####Cell Dimensions
# fig = plt.figure(figsize=fs_s_11)
# plt.plot(n_layers[1:],cell_a_diamond*1e10,zorder=2,color=rgbcode_diamond_m,linewidth=linewidth)
# plt.scatter(n_layers[1:],cell_a_diamond*1e10,zorder=3,color=rgbcode_black,edgecolor=rgbcode_diamond_m,linewidth=linewidth,marker=mark_d, s=marker_size)
# plt.hlines(y=ad0*1e10,xmin=n_layers[1],xmax=n_layers[-1],zorder=1,color=rgbcode_diamond_l,linewidth=linewidth,linestyle='dashed')
# plt.plot(n_layers[1:],cell_c_diamond*1e10,zorder=4,color=rgbcode_diamond_d,linewidth=linewidth)
# plt.scatter(n_layers[1:],cell_c_diamond*1e10,zorder=5,color=rgbcode_white,edgecolor=rgbcode_diamond_d,linewidth=linewidth,marker=mark_d, s=marker_size)
# plt.xlabel('N-ML',fontsize=axis_fontsize)
# plt.ylabel(r'Cell Size (\r{A})',fontsize=axis_fontsize)
# plt.xticks(fontsize=tick_fontsize)
# plt.yticks(fontsize=tick_fontsize)
#
# lines = [Line2D([0], [0], color=rgbcode_diamond_m,markerfacecolor=rgbcode_black,marker=mark_d,linewidth=universal_linewidth),
#          Line2D([0], [0], color=rgbcode_diamond_d,markerfacecolor=rgbcode_white, marker=mark_d, linewidth=universal_linewidth),
#          Line2D([0], [0], color=rgbcode_diamond_l,linestyle='dashed',linewidth=universal_linewidth),
# ]
# labels = ['N-ML Cell A','N-ML Cell C','Bulk Cell Sizes']
#
# plt.legend(lines,labels,loc='lower right')
# plot(showplot,'nml_cell0s_diamond')
#
#
# ####Bulk Mudulous vs ground state volume
# fig = plt.figure(figsize=fs_s_11)
# idxd = np.argsort(bulk_modulus_diamond)
# plt.plot(bulk_modulus_diamond[idxd]*1e-9,volumes0_diamond[idxd]*1e30,zorder=2,color=rgbcode_diamond_m,linewidth=linewidth)
# plt.scatter(bulk_modulus_diamond[idxd]*1e-9,volumes0_diamond[idxd]*1e30,zorder=3,color=rgbcode_diamond_m, linewidth=linewidth, s=marker_size,marker=mark_d)
# plt.hlines(y=vol0_d0*1e30,xmin=bulk_modulus_diamond[idxd][0]*1e-9,xmax=bulk_modulus_diamond[idxd][-1]*1e-9,zorder=1,color=adjust_lightness(rgbcode_diamond,adj_l),linewidth=linewidth,linestyle='dashed')
# plt.xlabel(r'${K}_{0}$ (GPa)',fontsize=axis_fontsize)
# plt.ylabel(r'${\mathit{V}}_{0}$ (\r{A}$^3$/atom)',fontsize=axis_fontsize)
# plt.xticks(fontsize=tick_fontsize)
# plt.yticks(fontsize=tick_fontsize)
#
# lines = [Line2D([0], [0], color=rgbcode_diamond_m,marker=mark_d,linewidth=universal_linewidth),
#          Line2D([0], [0], color=adjust_lightness(rgbcode_diamond,adj_l),linestyle='dashed',linewidth=universal_linewidth)]
# labels = ['N-ML','Bulk']
#
# plt.legend(lines,labels,loc='lower right')
# plot(showplot,'nml_bulkmod_volumes0_diamond')
#
#
# #####Bulk Mudulous Prime
# fig = plt.figure(figsize=fs_s_11)
# plt.plot(n_layers[1:],bulk_modulus_prime_diamond*1e-9,zorder=2,color=rgbcode_diamond_m,linewidth=linewidth)
# plt.scatter(n_layers[1:],bulk_modulus_prime_diamond*1e-9,zorder=3,color=rgbcode_diamond_m,linewidth=linewidth, s=marker_size, marker=mark_d)
# plt.hlines(y=bulk_modp0_d0*1e-9,xmin=n_layers[1],xmax=n_layers[-1],zorder=1,color=adjust_lightness(rgbcode_diamond,adj_l),linewidth=linewidth,linestyle='dashed')
# plt.xlabel('N-ML',fontsize=axis_fontsize)
# plt.ylabel(r'${{K}_{0}}\prime$',fontsize=axis_fontsize)
# plt.xticks(fontsize=tick_fontsize)
# plt.yticks(fontsize=tick_fontsize)
#
# lines = [Line2D([0], [0], color=rgbcode_diamond_m,marker=mark_d,linewidth=universal_linewidth),
#          Line2D([0], [0], color=adjust_lightness(rgbcode_diamond,adj_l),linestyle='dashed',linewidth=universal_linewidth)]
# labels = ['N-ML','Bulk']
#
# plt.legend(lines,labels,loc='upper right')
# plot(showplot,'nml_bulkmodulusprime_diamond')
#
#
#
#
# ####################################
# #######Betasn Plots
# ###Cell Dimensions
# fig = plt.figure(figsize=fs_s_11)
# plt.plot(n_layers[1:],cell_a_betasn*1e10,zorder=3,color=rgbcode_betasn_m,linewidth=linewidth)
# plt.scatter(n_layers[1:],cell_a_betasn*1e10,zorder=4,color=rgbcode_black,edgecolor=rgbcode_betasn_m,linewidth=linewidth,marker=mark_b, s=marker_size)
# plt.hlines(y=ab0*1e10,xmin=n_layers[1],xmax=n_layers[-1],zorder=1,color=adjust_lightness(rgbcode_betasn_m,adj_l),linewidth=linewidth,linestyle='dashed')
# plt.plot(n_layers[1:],cell_c_betasn*1e10,zorder=5,color=rgbcode_betasn_d,linewidth=linewidth)
# plt.scatter(n_layers[1:],cell_c_betasn*1e10,zorder=6,color=rgbcode_white,edgecolor=rgbcode_betasn_d,linewidth=linewidth,marker=mark_b, s=marker_size)
# plt.hlines(y=cb0*1e10,xmin=n_layers[1],xmax=n_layers[-1],zorder=2,color=adjust_lightness(rgbcode_betasn_d,adj_l),linewidth=linewidth,linestyle='dashed')
# plt.xlabel('N-ML',fontsize=axis_fontsize)
# plt.ylabel(r'Cell Size (\r{A})',fontsize=axis_fontsize)
# plt.xticks(fontsize=tick_fontsize)
# plt.yticks(fontsize=tick_fontsize)
#
#
# lines = [Line2D([0], [0], color=rgbcode_betasn_m,markerfacecolor=rgbcode_black,marker=mark_b,linewidth=universal_linewidth),
#          Line2D([0], [0], color=adjust_lightness(rgbcode_betasn_m,adj_l), linestyle='dashed', linewidth=universal_linewidth),
#          Line2D([0], [0], color=rgbcode_betasn_d,markerfacecolor=rgbcode_white, marker=mark_b, linewidth=universal_linewidth),
#          Line2D([0], [0], color=adjust_lightness(rgbcode_betasn_d,adj_l),linestyle='dashed',linewidth=universal_linewidth),
# ]
# labels = ['N-ML Cell A','Bulk Cell A','N-ML Cell C','Bulk Cell C']
#
# plt.legend(lines,labels,loc= 'lower right')
# plot(showplot,'nml_cell0s_betasn')
#
#
# ###Bulk Mudulous vs ground state volume
# fig = plt.figure(figsize=fs_s_11)
# idxb = np.argsort(bulk_modulus_betasn)
#
# plt.plot(bulk_modulus_betasn[idxb]*1e-9,volumes0_betasn[idxb]*1e30,zorder=2,color=rgbcode_betasn_m,linewidth=linewidth)
# plt.scatter(bulk_modulus_betasn[idxb]*1e-9,volumes0_betasn[idxb]*1e30,zorder=3,color=rgbcode_betasn_m, linewidth=linewidth, s=marker_size,marker=mark_b)
# plt.hlines(y=vol0_b0*1e30,xmin=bulk_modulus_betasn[idxb][0]*1e-9,xmax=bulk_modulus_betasn[idxb][-1]*1e-9,zorder=1,color=adjust_lightness(rgbcode_betasn,adj_l),linewidth=linewidth,linestyle='dashed')
# plt.xlabel(r'${K}_{0}$ (GPa)',fontsize=axis_fontsize)
# plt.ylabel(r'${V}_{0}$ (\r{A}$^3$/atom)',fontsize=axis_fontsize)
# plt.xticks(fontsize=tick_fontsize)
# plt.yticks(fontsize=tick_fontsize)
#
# lines = [Line2D([0], [0], color=rgbcode_betasn_m,marker=mark_b,linewidth=universal_linewidth),
#          Line2D([0], [0], color=adjust_lightness(rgbcode_betasn,adj_l),linestyle='dashed',linewidth=universal_linewidth)]
# labels = ['N-ML','Bulk']
#
# plt.legend(lines,labels,loc='lower left')
# plot(showplot,'nml_bulkmod_volumes0_betasn')
#
#
# ####Bulk Mudulous Prime
# fig = plt.figure(figsize=fs_s_11)
# plt.plot(n_layers[1:],bulk_modulus_prime_betasn*1e-9,zorder=2,color=rgbcode_betasn_m,linewidth=linewidth)
# plt.scatter(n_layers[1:],bulk_modulus_prime_betasn*1e-9,zorder=3,color=rgbcode_betasn_m,linewidth=linewidth, s=marker_size, marker=mark_b)
# plt.hlines(y=bulk_modp0_b0*1e-9,xmin=n_layers[1],xmax=n_layers[-1],zorder=1,color=adjust_lightness(rgbcode_betasn,adj_l),linewidth=linewidth,linestyle='dashed')
# plt.xlabel('N-ML',fontsize=axis_fontsize)
# plt.ylabel(r'${{K}_{0}}\prime$',fontsize=axis_fontsize)
# plt.xticks(fontsize=tick_fontsize)
# plt.yticks(fontsize=tick_fontsize)
#
# lines = [Line2D([0], [0], color=rgbcode_betasn_m,marker=mark_d,linewidth=universal_linewidth),
#          Line2D([0], [0], color=adjust_lightness(rgbcode_betasn,adj_l),linestyle='dashed',linewidth=universal_linewidth)]
# labels = ['N-ML','Bulk']
#
# plt.legend(lines,labels,loc='lower right')
# plot(showplot,'nml_bulkmodulusprime_betasn')
#
#
# ####Bulk Mudulous vs pressure
# fig = plt.figure(figsize=fs_s_11)
# plt.plot(-transition_pressures*1e-9,bulk_modulus_betasn*1e-9,zorder=2,color=rgbcode_betasn_m,linewidth=linewidth)
# plt.scatter(-transition_pressures*1e-9,bulk_modulus_betasn*1e-9,zorder=3,color=rgbcode_betasn_m,linewidth=linewidth, s=marker_size, marker=mark_b)
# plt.hlines(y=bulk_mod0_b0*1e-9,xmin=-transition_pressures[0]*1e-9,xmax=-transition_pressures[-1]*1e-9,zorder=1,color=adjust_lightness(rgbcode_betasn,adj_l),linewidth=linewidth,linestyle='dashed')
# plt.xlabel(r'${P}_{T}$ (GPa)',fontsize=axis_fontsize)
# plt.ylabel(r'${K}_{0}$ (GPa)',fontsize=axis_fontsize)
# plt.xticks(fontsize=tick_fontsize)
# plt.yticks(fontsize=tick_fontsize)
#
# lines = [Line2D([0], [0], color=rgbcode_betasn_m,marker=mark_d,linewidth=universal_linewidth),
#          Line2D([0], [0], color=adjust_lightness(rgbcode_betasn,adj_l),linestyle='dashed',linewidth=universal_linewidth)]
# labels = ['N-ML','Bulk']
#
# plt.legend(lines,labels,loc='lower right')
# plot(showplot,'nml_bulkmodulus_tpressure_betasn')
#
#
# ####Cohesive Energies
# fig = plt.figure(figsize=fs_s_11)
# plt.plot(n_layers[1:],cohesive_energy0_betasn*eV_per_joule,zorder=2,color=rgbcode_betasn_m,linewidth=linewidth)
# plt.scatter(n_layers[1:],cohesive_energy0_betasn*eV_per_joule,zorder=3,color=rgbcode_betasn_m,linewidth=linewidth, s=marker_size,marker=mark_b)
# plt.hlines(y=cohesive_energy0_b0*eV_per_joule,xmin=n_layers[1],xmax=n_layers[-1],zorder=1,color=adjust_lightness(rgbcode_betasn,adj_l),linewidth=linewidth,linestyle='dashed')
# plt.xlabel('N-ML',fontsize=axis_fontsize)
# plt.ylabel(r'$\mathit{E}$ (eV/atom)',fontsize=axis_fontsize)
# plt.xticks(fontsize=tick_fontsize)
# plt.yticks(fontsize=tick_fontsize)
#
# lines = [Line2D([0], [0], color=rgbcode_betasn_m,marker=mark_b,linewidth=universal_linewidth),
#          Line2D([0], [0], color=adjust_lightness(rgbcode_betasn,adj_l),linestyle='dashed',linewidth=universal_linewidth)]
# labels = ['N-ML','Bulk']
#
# plt.legend(lines,labels,loc='upper right')
# plot(showplot,'nml_cenergies_betasn')
# ####################################################
#
#
#
#
#
# fig = plt.figure(figsize=fs_s_11)
# plt.plot(bulk_modulus_diamond*1e-9,bonds_diamond,zorder=2,color=rgbcode_diamond_m,linewidth=linewidth)
# plt.scatter(bulk_modulus_diamond*1e-9,bonds_diamond,zorder=3,color=rgbcode_diamond_m, linewidth=linewidth, s=marker_size,marker=mark_b)
# plt.hlines(y=max_bond_d0,xmin=bulk_modulus_diamond[0]*1e-9,xmax=bulk_modulus_diamond[-1]*1e-9,zorder=1,color=adjust_lightness(rgbcode_diamond,adj_l),linewidth=linewidth,linestyle='dashed')
# plt.xlabel(r'${K}_{0}$ (GPa)',fontsize=axis_fontsize)
# plt.ylabel(r'${\mathit{V}}_{0}$ (\r{A}$^3$/atom)',fontsize=axis_fontsize)
# plt.xticks(fontsize=tick_fontsize)
# plt.yticks(fontsize=tick_fontsize)
#
# lines = [Line2D([0], [0], color=rgbcode_diamond_m,marker=mark_b,linewidth=universal_linewidth),
#          Line2D([0], [0], color=adjust_lightness(rgbcode_diamond,adj_l),linestyle='dashed',linewidth=universal_linewidth)]
# labels = ['N-ML','Bulk']
#
# plt.legend(lines,labels,loc='upper right')
# plot(showplot,'nml_bulkmod_volumes0_diamond')
#
# ####bond leghts
# plt.plot(n_layers[1:],bonds_diamond)
# plt.hlines(max_bond_d0,xmin=n_layers[1],xmax=n_layers[-1])
# plt.show()
#
# plt.plot(n_layers[1:],bonds_betasn_min)
# plt.hlines(min_bond_b0,xmin=n_layers[1],xmax=n_layers[-1])
# plt.show()
#
# plt.plot(n_layers[1:],bonds_betasn_max)
# plt.hlines(max_bond_b0,xmin=n_layers[1],xmax=n_layers[-1])
# plt.show()
#
# plt.plot(n_layers[1:],(bonds_betasn_max+bonds_betasn_min)/2)
# plt.hlines((max_bond_b0+min_bond_b0)/2,xmin=n_layers[1],xmax=n_layers[-1])
# plt.show()


l_b = 0

volumes0_diamond_exp = cell_a_diamond*cell_a_diamond*((cell_c_diamond*(n_layers[1:]-1+0.75))+(l_b*bonds_diamond*1e-10))/(8*n_layers[1:])
plt.plot(n_layers[1:],volumes0_diamond_exp*1e30,marker='s',color='blue')
plt.plot(n_layers[1:],volumes0_diamond*1e30,marker='s',alpha=0.3,color='blue')

plt.hlines(y=vol0_d0*1e30,xmin=n_layers[1],xmax=n_layers[-1])
plt.show()


cell_ct_diamond = 8*transition_volumes_diamond/(cell_a_diamond*cell_a_diamond)
volumes_diamond_exp = cell_a_diamond*cell_a_diamond*((cell_ct_diamond*(n_layers[1:]-1+0.75))+(l_b*bonds_diamond*1e-10))/(8*n_layers[1:])
plt.plot(n_layers[1:],volumes_diamond_exp*1e30,marker='s',color='blue')
plt.plot(n_layers[1:],transition_volumes_diamond*1e30,marker='s',alpha=0.3,color='blue')
plt.hlines(y=vol_d0*1e30,xmin=n_layers[1],xmax=n_layers[-1])
plt.show()


volumes0_betasn_exp = cell_a_betasn*cell_a_betasn*((cell_c_betasn*(n_layers[1:]-1+0.875))+(l_b*bonds_betasn_max*1e-10))/(8*n_layers[1:])
plt.plot(n_layers[1:],volumes0_betasn_exp*1e30,marker='s',color='red')
plt.plot(n_layers[1:],volumes0_betasn*1e30,marker='s',color='red',alpha=0.3)
plt.hlines(y=vol0_b0*1e30,xmin=n_layers[1],xmax=n_layers[-1],color='red')
plt.show()


cell_ct_betasn = 8*transition_volumes_betasn/(cell_a_betasn*cell_a_betasn)
volumes_betasn_exp = cell_a_betasn*cell_a_betasn*((cell_ct_betasn*(n_layers[1:]-1+0.875))+(l_b*bonds_betasn_max*1e-10))/(8*n_layers[1:])
plt.plot(n_layers[1:],volumes_betasn_exp*1e30,marker='s',color='red')
plt.plot(n_layers[1:],transition_volumes_betasn*1e30,marker='s',alpha=0.3,color='red')
plt.hlines(y=vol_b0*1e30,xmin=n_layers[1],xmax=n_layers[-1])
plt.show()