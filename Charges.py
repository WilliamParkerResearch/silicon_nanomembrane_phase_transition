import numpy as np
import matplotlib.pyplot as plt
from importlib import import_module
from ReferenceFiles.plot_formating import *
from matplotlib.lines import Line2D
from ReferenceFiles.plot_formating import rgbcode_betasn, rgbcode_diamond, universal_linewidth, marker_size, mark_d, mark_b

showplot = False
phase = 'diamond'
# n_layers = np.array([1,2,3,8])
phase = 'betasn'
n_layers = np.array([1,2,3])



if phase == 'diamond':
    rgbcode_electric = rgbcode_diamond
    mark_e = mark_d
    ext = 'D'
    miny, maxy = -0.055, 0.035
    phase_loc = 'lower center'
elif phase == 'betasn':
    rgbcode_electric = rgbcode_betasn
    mark_e = mark_b
    ext = 'B'
    miny, maxy = -0.025, 0.035
    phase_loc = 'upper center'

mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.sans-serif'] = "cmr10"
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['axes.formatter.useoffset'] = False

def charges_data(N_ML,exchange_correlation = 'PBE',phase='diamond'):
    import numpy as np
    from ReferenceFiles.plot_formating import rgbcode_betasn,rgbcode_diamond, universal_linewidth,marker_size,mark_d,mark_b
    from scipy.interpolate import UnivariateSpline
    from ReferenceFiles.FunctionDefinitions import spline_prep

    if phase == 'diamond':
        rgbcode_electric = rgbcode_diamond
        mark_e = mark_d
        ext='D'
    elif phase == 'betasn':
        rgbcode_electric = rgbcode_betasn
        mark_e = mark_b
        ext='B'
    directoryofdata = 'DataFolder' + '.' + exchange_correlation + '.' + 'charges' + '.' + 'Data_Charge_' + str(N_ML) + 'L'+ext
    data = import_module(directoryofdata)
    directoryofdata0 = 'DataFolder' + '.' + exchange_correlation + '.' + 'charges' + '.' + 'Data_Charge_' + '0' + 'L'+ext
    data0 = import_module(directoryofdata0)

    zpos = data.zpos
    charge = data.charge
    zpos_mod, charge_mod = spline_prep(zpos, charge)
    zpos_sim = np.linspace(zpos_mod[0],zpos_mod[-1],num=1000)
    charge_sim = UnivariateSpline(zpos_mod,charge_mod,k=5,s=0.005)

    zpos0 = data0.zpos
    charge0 = data0.charge
    cellc0 = data0.cella*data0.cell_c_a
    zpos_0_mod = np.array([])
    charge0_mod = np.array([])

    if N_ML>0:
        for i in np.arange(N_ML):
            z_add = i*cellc0
            zpos0_ext = zpos0+z_add
            zpos_0_mod = np.append(zpos_0_mod,zpos0_ext)
            charge0_mod = np.append(charge0_mod,charge0)


    charge0_mean = np.mean(charge0_mod)
    norm_charge_sim = (charge_sim(zpos_sim) - charge0_mean) / charge0_mean
    norm_charge = (charge - charge0_mean) / charge0_mean

    zpos_diff = zpos[-1]-zpos[0]
    norm_zpos_sim = zpos_sim/zpos_diff
    norm_zpos = zpos/zpos_diff
    return norm_zpos, norm_charge, norm_zpos_sim, norm_charge_sim

fig = plt.figure(figsize=fs_m_12)

k=1
lines = np.array([])
labels = []

for i in n_layers:
    chd = charges_data(i,phase=phase)
    n_entries = len(n_layers)
    minalpha = 2.3
    topalpha = 0.75
    pdiv = (1 / (n_entries))
    alpha = minalpha + k * pdiv * (topalpha - minalpha)
    rgbacode = adjust_lightness(rgbcode_electric, alpha)

    nline = Line2D([0], [0], color=rgbacode,marker=mark_e,linewidth=universal_linewidth)
    nlabel = str(i)+'-ML'

    lines = np.append(lines,nline)
    labels = np.append(labels,nlabel)

    if phase == 'diamond':
        if i == 1:
            zk = 3
            chd[1][0] = chd[1][0] + 0.0005
            chd[1][-1] = chd[1][-1] + 0.0005
        elif i == 2:
            zk = 1
            chd[1][0] = chd[1][0] + 0.00075
            chd[1][-1] = chd[1][-1] + 0.00075
        elif i == 3:
            zk = 2
        else:
            zk = k
            hgt=0.0


    if phase == 'betasn':
        if i == 1:
            zk = 2
            chd[1][0] = chd[1][0] + 0.0005
            chd[1][-1] = chd[1][-1] + 0.0005
        elif i == 2:
            zk = 1
        else:
            zk =k
            hgt=0.0



    plt.plot(chd[2],chd[3], color=rgbacode,linewidth=universal_linewidth,zorder=2*(zk-1)+1)
    plt.bar(chd[0],chd[1], color=rgbacode,width=0.03,zorder=2*(zk-1)+2,label = nlabel)
    # plt.scatter(chd[0],chd[1], color=rgbacode, s=marker_size, marker=mark_e,zorder=2*(k-1)+2)
    k+=1


s_pad = 0.05

plt.hlines(y=0,xmin=0-s_pad,xmax=1+s_pad,zorder=0,color=adjust_lightness(rgbcode_electric,0.5),linewidth=0.75*universal_linewidth,linestyle='dashed')
plt.xlabel(r'${Z}_{N-ML}/\Delta{Z}_{N-ML}$',fontsize=axis_fontsize)
plt.ylabel(r'$\Delta {Q}_{atom}/{Q}_{Bulk}$',fontsize=axis_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.ylim(miny,maxy)
plt.xlim(0-s_pad,1+s_pad)
plt.legend(loc=phase_loc)
plot(showplot,'charge_'+phase+'bar')