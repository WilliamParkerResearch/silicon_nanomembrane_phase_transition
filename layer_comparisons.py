import numpy as np
import matplotlib.pyplot as plt
from eos_information import eos_properties
from ReferenceFiles.cohesive_energy_values import *
from ReferenceFiles.unit_conversions import *
from scipy.optimize import curve_fit, fmin
from ReferenceFiles.FunctionDefinitions import mid, square_differences
from ReferenceFiles.plot_formating import *
from matplotlib.lines import Line2D
from ReferenceFiles.vasp_python_converter import vasp_data_modifier, bulk_cell_vstacker, ind_atom_distances, ztest
from ReferenceFiles.pair_distribution import format_pdf, plot_pdf
from scipy.signal import find_peaks
# <editor-fold desc="mpl mods">
#Use TeX fonts
# mpl.rcParams['text.usetex'] = True
# mpl.rcParams['font.sans-serif'] = "cmr10"
# mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['axes.formatter.useoffset'] = False
# </editor-fold>
# <editor-fold desc="plot format variables">
adj_l = 1.2
linewidth = universal_linewidth
# </editor-fold>


XC = 'PBE'
fold = 'eosztest'
prd = 1
prb = prd
showplot=True
if fold == 'eostest' or fold == 'eosztest':
    n_layers = np.array([0, 2, 4, 6, 8])
else:
    n_layers = np.arange(0, 9, 1)

############################################
############################################
# <editor-fold desc="Arrays for calculations">
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
bonds_diamond = np.array([])
bonds_betasn_min = np.array([])
bonds_betasn_max = np.array([])
bnddmin_arr = np.array([])
bndbmin_arr = np.array([])
bnddmax_arr = np.array([])
bndbmax_arr = np.array([])
# </editor-fold>

# <editor-fold desc="Property calculations">
for i in n_layers:
    minimum_distance=1.5
    maximum_distance = 6
#######################################################################################################################
    atomic_positions_d = vasp_data_modifier(i, 'diamond')[0]
    distance_values_d, pair_distribution_values_d = format_pdf(atomic_positions_d, mass_density=2300,
                                                                 minimum_distance=minimum_distance,
                                                                 maximum_distance=maximum_distance)
    br_idx_d = np.where(np.logical_and(distance_values_d >= 2.2, distance_values_d <= 2.5))[0]
    dv1d = distance_values_d[br_idx_d]
    pdv1d = pair_distribution_values_d[br_idx_d]
    maxbd_idx_d = np.argwhere(pdv1d == np.amax(pdv1d))[0][0]
    atom_bonds_mind = ind_atom_distances(1,atomic_positions_d)[0]
    atom_bonds_maxd = ind_atom_distances((i*8),atomic_positions_d)[0]
    cell_c_add = (atom_bonds_maxd+atom_bonds_mind)/2


    atomic_positions_b = vasp_data_modifier(i, 'betasn')[0]
    distance_values_b, pair_distribution_values_b = format_pdf(atomic_positions_b, 3370,
                                                                 minimum_distance=minimum_distance,
                                                                 maximum_distance=maximum_distance)
    br_idx_b = np.where(np.logical_and(distance_values_b >= 2.3, distance_values_d <= 2.85))[0]
    dv1b = distance_values_b[br_idx_b]
    pdv1b = pair_distribution_values_b[br_idx_b]
    maxbd_idx_b = np.argwhere(pdv1b == np.amax(pdv1b))[0][0]



    peaksd, _ = find_peaks(pdv1d, height=0)
    peaksb, _ = find_peaks(pdv1b, height=0,distance=150)
    if len(peaksd) > 1:
        pidxd = np.argsort(pdv1d[peaksd])
        peaksd = peaksd[pidxd[-1:]]

    if len(peaksb) > 2:
        pidxb = np.argsort(pdv1b[peaksb])
        peaksb = peaksb[pidxb[-2:]]

    max_bond_d = dv1d[peaksd][0]
    max_bond_b = np.sort(dv1b[peaksb])[-1]
    min_bond_b = np.sort(dv1b[peaksb])[0]
########################################################################################################################
    if i == 0:
        t_pres0, vol_d0, vol_b0, fitp_d0, fitp_b0,ad0,ab0,cd0,cb0,_,_,bnddmin0,bndbmin0,bnddmax0,bndbmax0 = eos_properties(i,exchange_correlation=XC,eos_folder=fold,propd=prd,propb=prb)
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

    else:
        t_pres, vol_d, vol_b, fitp_d, fitp_b,ad,ab,cd,cb,_,_,bnddmin,bndbmin,bnddmax,bndbmax= eos_properties(i,exchange_correlation=XC,eos_folder=fold,propd=prd,propb=prb)
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
        bnddmin_arr = np.append(bnddmin_arr,bnddmin)
        bndbmin_arr = np.append(bndbmin_arr,bndbmin)
        bnddmax_arr = np.append(bnddmax_arr,bnddmax)
        bndbmax_arr = np.append(bndbmax_arr,bndbmax)
# </editor-fold>

# <editor-fold desc="Cohesive Energy Conversions">
cohesive_energy0_diamond = -((energy0_diamond)-atom_energy_pbe)
cohesive_energy0_d0 = -((energy0_d0)-atom_energy_pbe)
cohesive_energy0_betasn = -((energy0_betasn)-atom_energy_pbe)
cohesive_energy0_b0 = -((energy0_b0)-atom_energy_pbe)
# </editor-fold>

# <editor-fold desc="Pressure smooth curve">
def property_curve_fitter(NML,x0,c1,c2):
    return x0+c1/(NML+c2)


def property_curve_fitter2(parameters,NML):
    return parameters[0]+parameters[1]/(NML+parameters[2])

optimized_parameters, covariance = curve_fit(property_curve_fitter, n_layers[1:],-transition_pressures,p0=[-t_pres0,t_pres0,0.01])

optimized_parameters2 = fmin(square_differences, [-t_pres0,t_pres0,0],
                                 args=(n_layers[1:], -transition_pressures, property_curve_fitter2), maxiter=100000)


m_layer = 8.25
cont_nml = np.linspace(0.75,m_layer,1000)
cont_transition_pressures = property_curve_fitter(cont_nml,optimized_parameters[0],optimized_parameters[1],optimized_parameters[2])
# </editor-fold>
############################################
############################################

########################################################################################################################
##### Plots ############################################################################################################
########################################################################################################################


####Cell C Dimensions Test
# # <editor-fold desc="Cell C Test">
# fig = plt.figure(figsize=fs_s_11)
# plt.plot(n_layers[1:],(cell_c_diamond*n_layers[1:])*1e10,zorder=4,color=rgbcode_diamond_d,linewidth=linewidth)
# plt.scatter(n_layers[1:],(cell_c_diamond*n_layers[1:])*1e10,zorder=5,color=rgbcode_white,edgecolor=rgbcode_diamond_d,linewidth=linewidth,marker=mark_d, s=marker_size)
# plt.plot(n_layers[1:],(cd0*(n_layers[1:]-1+0.75))*1e10,zorder=4,color=adjust_lightness(rgbcode_diamond_d,adj_l),linewidth=linewidth,linestyle='dashed')
#
#
# plt.plot(n_layers[1:],(cell_c_betasn*n_layers[1:])*1e10,zorder=4,color=rgbcode_betasn_d,linewidth=linewidth)
# plt.scatter(n_layers[1:],(cell_c_betasn*n_layers[1:])*1e10,zorder=5,color=rgbcode_white,edgecolor=rgbcode_betasn_d,linewidth=linewidth,marker=mark_d, s=marker_size)
# plt.plot(n_layers[1:],(cb0*(n_layers[1:]-1+0.875))*1e10,zorder=4,color=adjust_lightness(rgbcode_betasn_d,adj_l),linewidth=linewidth,linestyle='dashed')
#
#
# plt.xlabel('N-ML',fontsize=axis_fontsize)
# plt.ylabel(r'Cell Size (\r{A})',fontsize=axis_fontsize)
# plt.xticks(fontsize=tick_fontsize)
# plt.yticks(fontsize=tick_fontsize)
#
# plot(showplot,'nml_cell0s_test_diamond')
# # </editor-fold>


########################################################################################################################
##### Joint plots ######################################################################################################
########################################################################################################################


####Pressure Plot
# <editor-fold desc="P_T">
fig = plt.figure(figsize=fs_m_12)
# plt.plot(cont_nml,cont_transition_pressures*1e-9,zorder=2,color=adjust_lightness(rgbcode_pressure,0.9),linewidth=linewidth)
plt.plot(n_layers[1:],-transition_pressures*1e-9,zorder=3,color=adjust_lightness(rgbcode_pressure,0.9), linewidth=linewidth)
plt.scatter(n_layers[1:],-transition_pressures*1e-9,zorder=3,color=adjust_lightness(rgbcode_pressure,0.9), linewidth=linewidth, s=marker_size,marker=mark_p)
plt.hlines(y=-t_pres0*1e-9,zorder=4,xmin=n_layers[1],xmax=m_layer,color=adjust_lightness(rgbcode_pressure,adj_l),linewidth=linewidth,linestyle='dashed')
plt.hlines(y=0,zorder=1,xmin=n_layers[1],xmax=m_layer,color=rgbcode_black+(1,))


plt.xlabel('N-ML',fontsize=axis_fontsize)
plt.ylabel(r'$P_{T}$ (GPa)',fontsize=axis_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.yticks(fontsize=tick_fontsize)


lines = [Line2D([0], [0], color=adjust_lightness(rgbcode_pressure,0.9),marker=mark_p,linewidth=universal_linewidth),
         Line2D([0], [0], color=rgbcode_pressure,linestyle='dashed',linewidth=universal_linewidth)]
labels = ['N-ML','Bulk']


plt.legend(lines,labels,loc='lower right')
plot(showplot,'nml_tpressure')
# </editor-fold>

####Cohesive Energies
# <editor-fold desc="E_coh">
fig = plt.figure(figsize=fs_m_12)
plt.plot(n_layers[1:],cohesive_energy0_diamond*eV_per_joule,zorder=3,color=rgbcode_diamond_m,linewidth=linewidth)
plt.scatter(n_layers[1:],cohesive_energy0_diamond*eV_per_joule,zorder=4,color=rgbcode_diamond_m,linewidth=linewidth, s=marker_size,marker=mark_d)
plt.hlines(y=cohesive_energy0_d0*eV_per_joule,xmin=n_layers[1],xmax=n_layers[-1],zorder=2,color=adjust_lightness(rgbcode_diamond,adj_l),linewidth=linewidth,linestyle='dashed')

plt.plot(n_layers[1:],cohesive_energy0_betasn*eV_per_joule,zorder=5,color=rgbcode_betasn_m,linewidth=linewidth)
plt.scatter(n_layers[1:],cohesive_energy0_betasn*eV_per_joule,zorder=6,color=rgbcode_betasn_m,linewidth=linewidth, s=marker_size,marker=mark_b)
plt.hlines(y=cohesive_energy0_b0*eV_per_joule,xmin=n_layers[1],xmax=n_layers[-1],zorder=1,color=adjust_lightness(rgbcode_betasn,adj_l),linewidth=linewidth,linestyle='dashed')


plt.xlabel('N-ML',fontsize=axis_fontsize)
plt.ylabel(r'${E}_{coh}$ (eV/atom)',fontsize=axis_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.yticks(fontsize=tick_fontsize)

lines = [Line2D([0], [0], color=rgbcode_diamond_m,marker=mark_d,linewidth=universal_linewidth),
         Line2D([0], [0], color=rgbcode_betasn_m, marker=mark_b, linewidth=universal_linewidth),
         Line2D([0], [0], color=adjust_lightness(rgbcode_diamond,adj_l),linestyle='dashed',linewidth=universal_linewidth),
         Line2D([0], [0], color=adjust_lightness(rgbcode_betasn,adj_l),linestyle='dashed',linewidth=universal_linewidth)]
labels = ['N-ML Diamond*',r'N-ML $\beta$-tin*',' Bulk Diamond',r'Bulk $\beta$-tin']

plt.legend(lines,labels,loc='lower right')
plot(showplot,'nml_cenergies')
# </editor-fold>

########################################################################################################################
###### Diamond Properties ##############################################################################################
########################################################################################################################


####Ground State Volumes
# <editor-fold desc="Vol_0 Diamond">
fig = plt.figure(figsize=fs_s_11)
plt.plot(n_layers[1:],volumes0_diamond*1e30,zorder=2,color=adjust_lightness(rgbcode_diamond_m,0.9),linewidth=linewidth)
plt.scatter(n_layers[1:],volumes0_diamond*1e30,zorder=3,color=adjust_lightness(rgbcode_diamond_m,0.9), linewidth=linewidth, s=marker_size,marker=mark_d)
plt.hlines(y=vol0_d0*1e30,xmin=n_layers[1],xmax=n_layers[-1],zorder=1,color=adjust_lightness(rgbcode_diamond,adj_l-0.1),linewidth=linewidth,linestyle='dashed')
plt.xlabel('N-ML',fontsize=axis_fontsize)
plt.ylabel(r'${\mathit{V}}_{0}$ (\r{A}$^3$/atom)',fontsize=axis_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.ylim(19,25.5)

lines = [Line2D([0], [0], color=adjust_lightness(rgbcode_diamond_m,0.9),marker=mark_d,linewidth=universal_linewidth),
         Line2D([0], [0], color=adjust_lightness(rgbcode_diamond,adj_l-0.1),linestyle='dashed',linewidth=universal_linewidth)]
labels = ['N-ML','Bulk']

plt.legend(lines,labels,loc='upper right')
plot(showplot,'nml_volumes0_diamond')
# </editor-fold>

####Transition Volumes
# <editor-fold desc="Vol_T Diamond">
fig = plt.figure(figsize=fs_s_11)
plt.plot(n_layers[1:],transition_volumes_diamond*1e30,zorder=2,color=adjust_lightness(rgbcode_diamond_m,1.1),linewidth=linewidth)
plt.scatter(n_layers[1:],transition_volumes_diamond*1e30,zorder=3,color=adjust_lightness(rgbcode_diamond_m,1.1), linewidth=linewidth, s=marker_size,marker=mark_d)
plt.hlines(y=vol_d0*1e30,xmin=n_layers[1],xmax=n_layers[-1],zorder=1,color=adjust_lightness(rgbcode_diamond,adj_l+0.1),linewidth=linewidth,linestyle='dashed')
plt.xlabel('N-ML',fontsize=axis_fontsize)
plt.ylabel(r'${\mathit{V}}_{T}$ (\r{A}$^3$/atom)',fontsize=axis_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.ylim(19,25.5)

lines = [Line2D([0], [0], color=adjust_lightness(rgbcode_diamond_m,1.1),marker=mark_d,linewidth=universal_linewidth),
         Line2D([0], [0], color=adjust_lightness(rgbcode_diamond,adj_l+0.1),linestyle='dashed',linewidth=universal_linewidth)]
labels = ['N-ML','Bulk']

plt.legend(lines,labels,loc='upper right')
plot(showplot,'nml_tvolumes_diamond')
# </editor-fold>

####Cell Dimensions
# <editor-fold desc="CellSize Diamond">
fig = plt.figure(figsize=fs_s_11)
plt.plot(n_layers[1:],cell_a_diamond*1e10,zorder=2,color=rgbcode_diamond_m,linewidth=linewidth)
plt.scatter(n_layers[1:],cell_a_diamond*1e10,zorder=3,color=rgbcode_black,edgecolor=rgbcode_diamond_m,linewidth=linewidth,marker=mark_d, s=marker_size)
plt.hlines(y=ad0*1e10,xmin=n_layers[1],xmax=n_layers[-1],zorder=1,color=rgbcode_diamond_l,linewidth=linewidth,linestyle='dashed')

if fold == 'eos':
    plt.plot(n_layers[1:],(cell_c_diamond*1e10),zorder=5,color=rgbcode_diamond_d,linewidth=linewidth)
    plt.scatter(n_layers[1:],(cell_c_diamond*1e10),zorder=6,color=rgbcode_white,edgecolor=rgbcode_diamond_d,linewidth=linewidth,marker=mark_d, s=marker_size)
elif fold == 'eosz':
    plt.plot(n_layers[1:],(cell_c_diamond*1e10)/(n_layers[1:]),zorder=5,color=rgbcode_diamond_d,linewidth=linewidth)
    plt.scatter(n_layers[1:],(cell_c_diamond*1e10)/(n_layers[1:]),zorder=6,color=rgbcode_white,edgecolor=rgbcode_diamond_d,linewidth=linewidth,marker=mark_d, s=marker_size)
elif fold == 'eosztest' or fold == 'eostest':
    plt.plot(n_layers[1:],(cell_c_diamond*1e10)/(n_layers[1:]),zorder=5,color=rgbcode_diamond_d,linewidth=linewidth)
    plt.scatter(n_layers[1:],(cell_c_diamond*1e10)/(n_layers[1:]),zorder=6,color=rgbcode_white,edgecolor=rgbcode_diamond_d,linewidth=linewidth,marker=mark_d, s=marker_size)


# plt.plot(n_layers[1:],cell_c_diamond*1e10,zorder=4,color=rgbcode_diamond_d,linewidth=linewidth)
# plt.scatter(n_layers[1:],cell_c_diamond*1e10,zorder=5,color=rgbcode_white,edgecolor=rgbcode_diamond_d,linewidth=linewidth,marker=mark_d, s=marker_size)

plt.xlabel('N-ML',fontsize=axis_fontsize)
plt.ylabel(r'Cell Size (\r{A})',fontsize=axis_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.yticks(fontsize=tick_fontsize)

lines = [Line2D([0], [0], color=rgbcode_diamond_m,markerfacecolor=rgbcode_black,marker=mark_d,linewidth=universal_linewidth),
         Line2D([0], [0], color=rgbcode_diamond_d,markerfacecolor=rgbcode_white, marker=mark_d, linewidth=universal_linewidth),
         Line2D([0], [0], color=rgbcode_diamond_l,linestyle='dashed',linewidth=universal_linewidth),
]
labels = ['N-ML Cell A','N-ML Cell C/N-ML','Bulk Cell Sizes']

plt.legend(lines,labels,loc='upper right')
plot(showplot,'nml_cell0s_diamond')
# </editor-fold>

####Bulk Mudulous
# <editor-fold desc="K_0 Diamond">
fig = plt.figure(figsize=fs_s_11)
plt.plot(n_layers[1:],bulk_modulus_diamond*1e-9,zorder=2,color=rgbcode_diamond_m,linewidth=linewidth)
plt.scatter(n_layers[1:],bulk_modulus_diamond*1e-9,zorder=3,color=rgbcode_diamond_m,linewidth=linewidth, s=marker_size, marker=mark_d)
plt.hlines(y=bulk_mod0_d0*1e-9,xmin=n_layers[1],xmax=n_layers[-1],zorder=1,color=adjust_lightness(rgbcode_diamond,adj_l),linewidth=linewidth,linestyle='dashed')
plt.xlabel('N-ML',fontsize=axis_fontsize)
plt.ylabel(r'${K}_{0}$ (GPa)',fontsize=axis_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.ylim(140,365)


lines = [Line2D([0], [0], color=rgbcode_diamond_m,marker=mark_d,linewidth=universal_linewidth),
         Line2D([0], [0], color=adjust_lightness(rgbcode_diamond,adj_l),linestyle='dashed',linewidth=universal_linewidth)]
labels = ['N-ML','Bulk']

plt.legend(lines,labels,loc='upper right')
plot(showplot,'nml_bulkmodulus_diamond')
# </editor-fold>


########################################################################################################################
####### BetaSn Properties ##############################################################################################
########################################################################################################################


####Ground State Volumes
# <editor-fold desc="Vol_0 Betasn">
fig = plt.figure(figsize=fs_s_11)
plt.plot(n_layers[1:],volumes0_betasn*1e30,zorder=2,color=adjust_lightness(rgbcode_betasn_m,0.9),linewidth=linewidth)
plt.scatter(n_layers[1:],volumes0_betasn*1e30,zorder=3,color=adjust_lightness(rgbcode_betasn_m,0.9), linewidth=linewidth, s=marker_size,marker=mark_b)
plt.hlines(y=vol0_b0*1e30,xmin=n_layers[1],xmax=n_layers[-1],zorder=1,color=adjust_lightness(rgbcode_betasn,adj_l-0.1),linewidth=linewidth,linestyle='dashed')
plt.xlabel('N-ML',fontsize=axis_fontsize)
plt.ylabel(r'${\mathit{V}}_{0}$ (\r{A}$^3$/atom)',fontsize=axis_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.ylim(14.5,19)


lines = [Line2D([0], [0], color=adjust_lightness(rgbcode_betasn_m,0.9),marker=mark_b,linewidth=universal_linewidth),
         Line2D([0], [0], color=adjust_lightness(rgbcode_betasn,adj_l-0.1),linestyle='dashed',linewidth=universal_linewidth)]
labels = ['N-ML','Bulk']

plt.legend(lines,labels,loc='upper right')
plot(showplot,'nml_volumes0_betasn')
# </editor-fold>

####Transition Volumes
# <editor-fold desc="Vol_T Betasn">
fig = plt.figure(figsize=fs_s_11)
plt.plot(n_layers[1:],transition_volumes_betasn*1e30,zorder=2,color=adjust_lightness(rgbcode_betasn_m,1.1),linewidth=linewidth)
plt.scatter(n_layers[1:],transition_volumes_betasn*1e30,zorder=3,color=adjust_lightness(rgbcode_betasn_m,1.1), linewidth=linewidth, s=marker_size,marker=mark_b)
plt.hlines(y=vol_b0*1e30,xmin=n_layers[1],xmax=n_layers[-1],zorder=1,color=adjust_lightness(rgbcode_betasn,adj_l+0.1),linewidth=linewidth,linestyle='dashed')
plt.xlabel('N-ML',fontsize=axis_fontsize)
plt.ylabel(r'${\mathit{V}}_{T}$ (\r{A}$^3$/atom)',fontsize=axis_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.ylim(14.5,19)

lines = [Line2D([0], [0], color=adjust_lightness(rgbcode_betasn_m,1.1),marker=mark_b,linewidth=universal_linewidth),
         Line2D([0], [0], color=adjust_lightness(rgbcode_betasn,adj_l+0.1),linestyle='dashed',linewidth=universal_linewidth)]
labels = ['N-ML','Bulk']

plt.legend(lines,labels,loc='upper right')
plot(showplot,'nml_tvolumes_betasn')
# </editor-fold>

####Cell dimensions
# <editor-fold desc="CellSize Betasn">
fig = plt.figure(figsize=fs_s_11)
plt.plot(n_layers[1:],cell_a_betasn*1e10,zorder=3,color=rgbcode_betasn_m,linewidth=linewidth)
plt.scatter(n_layers[1:],cell_a_betasn*1e10,zorder=4,color=rgbcode_black,edgecolor=rgbcode_betasn_m,linewidth=linewidth,marker=mark_b, s=marker_size)
plt.hlines(y=ab0*1e10,xmin=n_layers[1],xmax=n_layers[-1],zorder=1,color=adjust_lightness(rgbcode_betasn_m,adj_l),linewidth=linewidth,linestyle='dashed')

if fold == 'eos':
    plt.plot(n_layers[1:],(cell_c_betasn*1e10/2),zorder=5,color=rgbcode_betasn_d,linewidth=linewidth)
    plt.scatter(n_layers[1:],(cell_c_betasn*1e10/2),zorder=6,color=rgbcode_white,edgecolor=rgbcode_betasn_d,linewidth=linewidth,marker=mark_b, s=marker_size)
elif fold == 'eosz':
    plt.plot(n_layers[1:],(cell_c_betasn*1e10)/(2*n_layers[1:]),zorder=5,color=rgbcode_betasn_d,linewidth=linewidth)
    plt.scatter(n_layers[1:],(cell_c_betasn*1e10)/(2*n_layers[1:]),zorder=6,color=rgbcode_white,edgecolor=rgbcode_betasn_d,linewidth=linewidth,marker=mark_b, s=marker_size)
elif fold == 'eosztest' or fold == 'eostest':
    plt.plot(n_layers[1:],(cell_c_betasn*1e10)/(n_layers[1:]),zorder=5,color=rgbcode_betasn_d,linewidth=linewidth)
    plt.scatter(n_layers[1:],(cell_c_betasn*1e10)/(n_layers[1:]),zorder=6,color=rgbcode_white,edgecolor=rgbcode_betasn_d,linewidth=linewidth,marker=mark_b, s=marker_size)

# plt.plot(n_layers[1:],(cell_c_betasn*1e10/2),zorder=5,color=rgbcode_betasn_d,linewidth=linewidth)
# plt.scatter(n_layers[1:],(cell_c_betasn*1e10/2),zorder=6,color=rgbcode_white,edgecolor=rgbcode_betasn_d,linewidth=linewidth,marker=mark_b, s=marker_size)

cellc0 = ((cb0*1e10/2))#+(bonds_betasn_max+bonds_betasn_min)/4)
plt.hlines(y=cellc0,xmin=n_layers[1],xmax=n_layers[-1],zorder=2,color=adjust_lightness(rgbcode_betasn_d,adj_l),linewidth=linewidth,linestyle='dashed')

plt.xlabel('N-ML',fontsize=axis_fontsize)
plt.ylabel(r'Cell Size (\r{A})',fontsize=axis_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.yticks(fontsize=tick_fontsize)


lines = [Line2D([0], [0], color=rgbcode_betasn_m,markerfacecolor=rgbcode_black,marker=mark_b,linewidth=universal_linewidth),
         Line2D([0], [0], color=adjust_lightness(rgbcode_betasn_m,adj_l), linestyle='dashed', linewidth=universal_linewidth),
         Line2D([0], [0], color=rgbcode_betasn_d,markerfacecolor=rgbcode_white, marker=mark_b, linewidth=universal_linewidth),
         Line2D([0], [0], color=adjust_lightness(rgbcode_betasn_d,adj_l),linestyle='dashed',linewidth=universal_linewidth),
]
labels = ['N-ML Cell A','Bulk Cell A','N-ML Cell C/N-ML','Bulk Cell C']

plt.legend(lines,labels) #bbox_to_anchor=(0.5, 0.4)
plot(showplot,'nml_cell0s_betasn')
# </editor-fold>

####Bulk Mudulous
# <editor-fold desc="K_0 Betasn">
fig = plt.figure(figsize=fs_s_11)
plt.plot(n_layers[1:],bulk_modulus_betasn*1e-9,zorder=2,color=rgbcode_betasn_m,linewidth=linewidth)
plt.scatter(n_layers[1:],bulk_modulus_betasn*1e-9,zorder=3,color=rgbcode_betasn_m,linewidth=linewidth, s=marker_size,marker=mark_b)
plt.hlines(y=bulk_mod0_b0*1e-9,xmin=n_layers[1],xmax=n_layers[-1],zorder=1,color=adjust_lightness(rgbcode_betasn,adj_l),linewidth=linewidth,linestyle='dashed')
plt.xlabel('N-ML',fontsize=axis_fontsize)
plt.ylabel(r'${K}_{0}$ (GPa)',fontsize=axis_fontsize)
plt.ylim(140,365)
plt.xticks(fontsize=tick_fontsize)
plt.yticks(fontsize=tick_fontsize)

lines = [Line2D([0], [0], color=rgbcode_betasn_m,marker=mark_b,linewidth=universal_linewidth),
         Line2D([0], [0], color=adjust_lightness(rgbcode_betasn,adj_l),linestyle='dashed',linewidth=universal_linewidth)]
labels = ['N-ML','Bulk']

plt.legend(lines,labels,loc='upper right')
plot(showplot,'nml_bulkmodulus_betasn')
# </editor-fold>