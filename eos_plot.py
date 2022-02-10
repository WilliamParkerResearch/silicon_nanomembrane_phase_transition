import numpy as np
import matplotlib.pyplot as plt
from ReferenceFiles.plot_formating import *
import matplotlib as mpl
from matplotlib.lines import Line2D
from eos_information import eos_properties

# Use TeX fonts
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.sans-serif'] = "cmr10"
mpl.rcParams['font.family'] = 'serif'

def plot_eos(N_ML,alpha=1,dw=1,ds=0,exchange_correlation='PBE', OF=False, etype='ce', n_etype=0, eos_type = 'vinet',zero=0,fold='eos'):
    import numpy as np
    import matplotlib.pyplot as plt
    from eos_information import eos_properties
    from ReferenceFiles.equations_of_state import birch_murnaghan, murnaghan, vinet, \
        pressure_from_energy_equation_of_state
    from ReferenceFiles.plot_formating import rgbcode_diamond, rgbcode_betasn,rgbcode_pressure, universal_linewidth, mark_d, mark_b,marker_size, adjust_lightness
    from importlib import import_module

    rgbcode_diamond = adjust_lightness(rgbcode_diamond,0.85)
    rgbcode_betasn = adjust_lightness(rgbcode_betasn,1)
    rgbcode_pressure = adjust_lightness(rgbcode_pressure,1.15)

    if OF == False:
        directoryofdata = 'DataFolder' + '.' + exchange_correlation + '.' + fold + '.' + 'Data_' + str(N_ML) + 'L'
    else:
        directoryofdata = 'DataFolder' + '.' + exchange_correlation + '.' + fold + '.' + 'Data_' + str(N_ML) + 'L' + '_' + etype + str(n_etype)

    data = import_module(directoryofdata)

    energies_d = data.total_energies_strain_diamond
    # volumes_d = data.volumes_sim_diamond
    energies_b = data.total_energies_strain_betasn
    # volumes_b = data.volumes_sim_betasn
    
    if eos_type == 'birch-murnaghan':
        def eos(p, v):
            return birch_murnaghan(p, v)
    elif eos_type == 'murnaghan':
        def eos(p, v):
            return murnaghan(p, v)
    elif eos_type == 'vinet':
        def eos(p, v):
            return vinet(p, v)

    rgbacode_diamond = rgbcode_diamond + (alpha,)
    rgbacode_betasn = rgbcode_betasn + (alpha,)
    rgbacode_pressure = rgbcode_pressure +(alpha,)
    linewidth = universal_linewidth
    per_reach = -0.05
    num_vol = 10000
    t_pres0, vol_d0, vol_b0, fitp_d0, fitp_b0,_,_,_,_,volumes_d0,volumes_b0 = eos_properties(1, eos_type=eos_type,eos_folder=fold,propd=1,propb=1)
    t_press, vol_d, vol_b, fitp_d, fitp_b,_,_,_,_,volumes_d,volumes_b = eos_properties(N_ML=N_ML, eos_type=eos_type,eos_folder=fold,propd=1,propb=1)

    delta_vol_0 = fitp_d0[-1]-fitp_b0[-1]

    prop = 1+per_reach

    vol_0_d0 = fitp_d0[-1]
    vol_0_b0 = fitp_b0[-1]
    vol_0_d = fitp_d[-1]
    vol_0_b = fitp_b[-1]


    min_vol_d = vol_0_d - prop*(delta_vol_0/2)
    max_vol_d = vol_0_d + prop*(delta_vol_0/2)
    min_vol_b = vol_0_b - prop*(delta_vol_0/2)
    max_vol_b = vol_0_b + prop*(delta_vol_0/2)

    E_i_d = eos(fitp_d,min_vol_d)
    E_f_d = eos(fitp_d,max_vol_d)
    E_i_b = eos(fitp_b,min_vol_b)
    E_f_b = eos(fitp_b,max_vol_b)
    E_min_edge_d = np.amin(np.array([E_i_d,E_f_d]))
    E_min_edge_b = np.amin(np.array([E_i_b,E_f_b]))


    min_vol_d0 = vol_0_d0 - prop*(delta_vol_0/2)
    max_vol_d0 = vol_0_d0 + prop*(delta_vol_0/2)
    min_vol_b0 = vol_0_b0 - prop*(delta_vol_0/2)
    max_vol_b0 = vol_0_b0 + prop*(delta_vol_0/2)

    E_i_d0 = eos(fitp_d0,min_vol_d0)
    E_f_d0 = eos(fitp_d0,max_vol_d0)
    E_i_b0 = eos(fitp_b0,min_vol_b0)
    E_f_b0 = eos(fitp_b0,max_vol_b0)
    E_min_edge_d0 = np.amin(np.array([E_i_d0,E_f_d0]))
    E_min_edge_b0 = np.amin(np.array([E_i_b0,E_f_b0]))

    volumes_diamond0 = np.linspace(min_vol_d0,max_vol_d0,num_vol)
    volumes_betasn0 = np.linspace(min_vol_b0,max_vol_b0,num_vol)
    eos_diamond0 = eos(fitp_d0,volumes_diamond0)
    eos_betasn0 = eos(fitp_b0,volumes_betasn0)

    volumes_diamond = np.linspace(min_vol_d,max_vol_d,num_vol)
    volumes_betasn = np.linspace(min_vol_b,max_vol_b,num_vol)
    eos_diamond = eos(fitp_d,volumes_diamond)
    eos_betasn = eos(fitp_b,volumes_betasn)

    idx_b0 = np.argwhere(eos_betasn0<=(E_min_edge_b0))
    volumes_betasn0 = volumes_betasn0[idx_b0]
    idx_d0 = np.argwhere(eos_diamond0<=(E_min_edge_d0))
    volumes_diamond0 = volumes_diamond0[idx_d0]


    idx_d = np.argwhere(eos_diamond<=(E_min_edge_d))
    volumes_diamond = volumes_diamond[idx_d]
    eos_diamond = eos_diamond[idx_d]
    idx_b = np.argwhere(eos_betasn<=(E_min_edge_b))
    volumes_betasn = volumes_betasn[idx_b]
    eos_betasn = eos_betasn[idx_b]

    idx_s_d = np.argwhere(energies_d<=(E_min_edge_d))
    energies_s_d = energies_d[idx_s_d]
    volumes_s_d = volumes_d[idx_s_d]
    idx_s_b = np.argwhere(energies_b<=(E_min_edge_b))
    energies_s_b = energies_b[idx_s_b]
    volumes_s_b = volumes_b[idx_s_b]
    


    volumes_col = np.append(volumes_betasn0,volumes_diamond0)
    # volumes = np.linspace(np.amin(volumes_col),np.amax(volumes_col))
    volumes = np.linspace(vol_b,vol_d)

    tangent = t_press*(volumes-vol_d) + eos(fitp_d,vol_d)
    if N_ML == 0:
        mark_color_d = 'white'
        edge_color_d = rgbacode_diamond
        mark_color_b = 'white'
        edge_color_b = rgbacode_betasn
    else:
        mark_color_d = rgbacode_diamond
        edge_color_d = rgbacode_diamond
        mark_color_b = rgbacode_betasn
        edge_color_b = rgbacode_betasn

    plt.plot(volumes*1e30,tangent/1.6e-19-zero, zorder=1 ,color=rgbacode_pressure, linewidth=linewidth, dashes=(dw,ds))
    plt.plot(volumes_diamond*1e30,eos_diamond/1.6e-19-zero, zorder=4, color=rgbacode_diamond, linewidth=linewidth, dashes=(dw,ds))
    plt.scatter(volumes_s_d*1e30,energies_s_d/1.6e-19-zero, zorder=5, color=mark_color_d, edgecolor=edge_color_d, linewidth=linewidth,marker=mark_d,s=marker_size)
    plt.plot(volumes_betasn*1e30,eos_betasn/1.6e-19-zero, zorder=2, color= rgbacode_betasn, linewidth=linewidth, dashes=(dw,ds))
    plt.scatter(volumes_s_b*1e30,energies_s_b/1.6e-19-zero, zorder=3, color=mark_color_b, edgecolor=edge_color_b, linewidth=linewidth, marker=mark_b,s=marker_size)
    return

alpha0=0.75
dw0,ds0 = 3,3

E_shift = eos_properties(0)[3][0]/1.6e-19
fig = plt.figure(figsize=fs_l_12)
plot_eos(0,alpha=alpha0,dw=dw0,ds=ds0,zero=E_shift,fold='eosz')
plot_eos(1,zero=E_shift,fold='eosz')
plt.xlabel(r'$V$ (\r{A}$^3$/atom)',fontsize=axis_fontsize)
plt.ylabel(r'$\mathit{E}$ (eV/atom)',fontsize=axis_fontsize)
plt.tick_params(axis='x', which='major', labelsize=tick_fontsize)
plt.tick_params(axis='y', which='major', labelsize=tick_fontsize)

lines = [Line2D([0], [0], color=rgbcode_diamond,marker=mark_d,linewidth=universal_linewidth),
         Line2D([0], [0], color=rgbcode_betasn, marker=mark_b, linewidth=universal_linewidth),
         Line2D([0], [0], color=rgbcode_diamond+(alpha0,),dashes=(dw0,ds0),marker=mark_d,markerfacecolor='w',linewidth=universal_linewidth),
         Line2D([0], [0], color=rgbcode_betasn+(alpha0,), dashes=(dw0, ds0), marker=mark_b, markerfacecolor='w',
                linewidth=universal_linewidth)
         ]
labels = ['1-ML Diamond*',r'1-ML $\beta$-tin*',' Bulk Diamond',r'Bulk $\beta$-tin']

plt.legend(lines,labels,loc='center right')
plot(False,'EOS')

