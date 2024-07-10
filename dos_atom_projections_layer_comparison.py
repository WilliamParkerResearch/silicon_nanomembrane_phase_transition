import matplotlib.pyplot as plt
from matplotlib import rcParamsDefault
from ReferenceFiles.plot_formating import *
import numpy as np
from matplotlib.lines import Line2D

phase = 'diamond'

if phase == 'diamond':
    rgbcode_electric = rgbcode_diamond
    mark_e = mark_d
elif phase == 'betasn':
    rgbcode_electric = rgbcode_betasn
    mark_e = mark_b


def dos_atom_proj_arrays(N_ML,phase = 'diamond',exchange_correlation = 'PBE'):
    import numpy as np
    from ReferenceFiles.dos_functions import pdos
    from ReferenceFiles.FunctionDefinitions import zpos_sort_idx, find_nearest,spline_prep
    from ReferenceFiles.vasp_python_converter import vasp_data_modifier
    from importlib import import_module



    N_ML = str(N_ML)


    order_idx = zpos_sort_idx(N_ML,phase)
    zpos = vasp_data_modifier(N_ML,phase)[0][:,2][order_idx]-vasp_data_modifier(N_ML,phase)[0][:,2][order_idx][0]
    nat = len(zpos)

    if phase == 'diamond':
        prefix = 'Si.Fd-3m_'
    elif phase == 'betasn':
        prefix = 'Si.I4_1amd_'


    directoryofdata='DataFolder'+'.'+exchange_correlation+'.'+'bands'+'.'+phase+'.'+'Data_Bands_'+str(N_ML)+'L'
    data = import_module(directoryofdata)

    fermi_energy = data.fermi_energy
    pref = 'DataFolder/'+exchange_correlation+'/pdos/'+phase + '/'+str(N_ML)+'ML/'+prefix+str(N_ML)+'ML.k.pdos_'




    et, pdost = pdos(pref+'tot','p')

    idx_0 = find_nearest(et-fermi_energy,0)[1]
    pdost_0 = pdost[idx_0]


    def atom_pdos(atoms):
        es1, pdoss1 = np.zeros([1,len(et)])[0], np.zeros([1,len(pdost)])[0]
        es2, pdoss2 = np.zeros([1,len(et)])[0], np.zeros([1,len(pdost)])[0]

        for i in atoms:
            atom_number = order_idx[i-1]+1
            e1, pdos1 = pdos(pref+'atm#'+str(atom_number)+'(Si)_wfc#1(s)')
            e2, pdos2 = pdos(pref+'atm#'+str(atom_number)+'(Si)_wfc#2(p)')
            pdoss1 = pdoss1 + pdos1
            pdoss2 = pdoss2 + pdos2
        atoms_pdoss = pdoss1+pdoss2
        return atoms_pdoss


    dos_0_arr = np.zeros(nat)

    for i in (np.arange(nat)+1):
        print(i)
        atoms_dos = atom_pdos(np.array([i]))[idx_0]
        dos_0_arr[i-1] = atoms_dos

    dos_0_per_arr = dos_0_arr/pdost_0



    from scipy.interpolate import UnivariateSpline,InterpolatedUnivariateSpline,LSQUnivariateSpline

    zpos_sim = np.linspace(zpos[0],zpos[-1],num=10000,endpoint=True)

    zpos_smod,dos_0_sarr = spline_prep(zpos,dos_0_per_arr)


    midzp = (zpos_smod[1:]+zpos_smod[:-1])/2


    k_val = 2

    if N_ML == '0':
        k_val = 2
    elif N_ML == '1':
        k_val = 2
        t = [(zpos_smod[0]+zpos_smod[-1])/2]
    elif N_ML == '2':
        t = [zpos_smod[2],zpos_smod[-2]]
    elif N_ML == '3':
        t = [zpos_smod[2],zpos_smod[-2]]
    elif N_ML == '8':
        t = [zpos_smod[5],zpos_smod[-2]]


    dos_0_per_arr_sim = LSQUnivariateSpline(zpos_smod,dos_0_sarr,t,k=2)

    zpos_diff = np.sort(zpos)[-1] - np.sort(zpos)[0]
    return zpos/zpos_diff, dos_0_per_arr, zpos_sim/zpos_diff, dos_0_per_arr_sim(zpos_sim)



n_layers = np.array([1,2,3,8])



fig = plt.figure(figsize=fs_m_12)

k=1
lines = np.array([])
labels = []

for i in n_layers:
    n_entries = len(n_layers)
    minalpha = 2.3
    topalpha = 0.75
    pdiv = (1 / (n_entries))
    alpha = minalpha + k * pdiv * (topalpha - minalpha)
    rgbacode = adjust_lightness(rgbcode_electric, alpha)

    zpos, dos_0_per_arr, zpos_sim, dos_0_per_arr_sim = dos_atom_proj_arrays(i)

    nline = Line2D([0], [0], color=rgbacode,marker=mark_e,linewidth=universal_linewidth)
    nlabel = str(i)+'-ML'

    lines = np.append(lines,nline)
    labels = np.append(labels,nlabel)


    if phase == 'diamond':
        if i == 1:
            zk = 1
        elif i == 2:
            zk = 2
            zpos[0] = zpos[0]-0.0005
            zpos[-1] = zpos[-1]+0.0005
        elif i == 3:
            zk = 2
        elif i == 8:
            zk = 4
            dos_0_per_arr[1] = dos_0_per_arr[0]+0.0005
            dos_0_per_arr[-1] = dos_0_per_arr+0.0005
    plt.plot(zpos_sim, dos_0_per_arr_sim, color=rgbacode, linewidth=universal_linewidth,zorder=2*(zk-1)+1)
    plt.bar(zpos, dos_0_per_arr,color=rgbacode, width=0.03,zorder=2*(zk-1)+2,label = nlabel)
    k+=1










plt.ylim(0,0.31)
plt.gca().set_ylim(bottom=0)
plt.xlabel(r'${Z}_{N-ML}/\Delta{Z}_{N-ML}$',fontsize=axis_fontsize)
plt.ylabel(r'${\rho}_{atom}/{\rho}_{total}  $',fontsize=axis_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.legend(bbox_to_anchor=(0.85, 0.55))
plot(True,name='pdos_atom_contributios_'+phase+'.png')
