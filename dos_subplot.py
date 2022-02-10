import matplotlib.pyplot as plt
from matplotlib import rcParamsDefault
import numpy as np
from ReferenceFiles.dos_functions import *
from ReferenceFiles.plot_formating import *
from ReferenceFiles.FunctionDefinitions import zpos_sort_idx
from importlib import import_module
import matplotlib.patches as mpatches

def dos_plotter(N_ML,exchange_correlation = 'PBE',phase = 'diamond',lightness = 1):
    N_ML = str(N_ML)

    order_idx = zpos_sort_idx(N_ML,phase)


    if N_ML == '0':
        nat = 8
    else:
        nat = int(float(N_ML)*8)

    if phase == 'diamond':
        prefix = 'Si.Fd-3m_'
        rgbcode_electric = rgbcode_diamond
        mark_e = mark_d
    elif phase == 'betasn':
        prefix = 'Si.I4_1amd_'
        rgbcode_electric = rgbcode_betasn
        mark_e = mark_b

    directoryofdata='DataFolder'+'.'+exchange_correlation+'.'+'bands'+'.'+phase+'.'+'Data_Bands_'+N_ML+'L'
    data = import_module(directoryofdata)
    fermi_energy = data.fermi_energy
    pref = 'DataFolder/'+exchange_correlation+'/pdos/'+phase + '/'+str(N_ML)+'ML/'+prefix+str(N_ML)+'ML.k.pdos_'




    et, pdost = pdos(pref+'tot','p')

    def atom_pdos(atoms):
        es1, pdoss1 = np.zeros([1,len(et)])[0], np.zeros([1,len(pdost)])[0]
        es2, pdoss2 = np.zeros([1,len(et)])[0], np.zeros([1,len(pdost)])[0]

        for i in atoms:
            print(str(i))
            atom_number = order_idx[i-1]+1
            e1, pdos1 = pdos(pref+'atm#'+str(atom_number)+'(Si)_wfc#1(s)')
            e2, pdos2 = pdos(pref+'atm#'+str(atom_number)+'(Si)_wfc#2(p)')
            pdoss1 = pdoss1 + pdos1
            pdoss2 = pdoss2 + pdos2
        atoms_pdoss = pdoss1+pdoss2
        return atoms_pdoss


    atoms_dos = atom_pdos((np.array([])))
    nel = nat*4

    if N_ML=='0':
        labelc = 'Bulk'
    else:
        labelc = str(N_ML) + '-ML'
    plt.fill_between(et - fermi_energy,pdost/nat,color=adjust_lightness(rgbcode_electric,lightness))
    plt.plot(et - fermi_energy,pdost/nat,color=adjust_lightness(rgbcode_electric,lightness),linewidth=universal_linewidth,label=labelc)

    plt.fill_between(et - fermi_energy,atoms_dos,color=adjust_lightness(rgbcode_electric,lightness-0.25))
    plt.plot(et - fermi_energy,atoms_dos,color=adjust_lightness(rgbcode_electric,lightness-0.4),linewidth=universal_linewidth)

    plt.hlines(y=0,xmin=np.amin(et)-fermi_energy,xmax=np.amax(et)-fermi_energy,color=rgbcode_black,linewidth=universal_linewidth)
    plt.xlim(-5,5)
    plt.gca().set_ylim(bottom=0)
    plt.xticks(fontsize=tick_fontsize)
    plt.yticks(fontsize=tick_fontsize)
    plt.yticks([])
    return


phase = 'diamond'


fig = plt.figure(figsize=fs_s_23)
fig.subplots_adjust(hspace=0)


ax1 = plt.subplot(212)
dos_plotter(1,phase=phase,lightness=1.15)
color1 = plt.gca().lines[0].get_color()
label1 = plt.gca().lines[0].get_label()
plt.xlabel(r'$\varepsilon_{\rm KS}$ (eV)', fontsize=axis_fontsize)


ax0 = plt.subplot(211)
dos_plotter(0,phase=phase,lightness=0.85)
color0 = plt.gca().lines[0].get_color()
label0 = plt.gca().lines[0].get_label()
ax0.axes.xaxis.set_visible(False)


DOS0 = mpatches.Patch(color=color0, label=label0)
DOS1 = mpatches.Patch(color=color1, label=label1)

plt.legend(handles=[DOS0,DOS1], loc='upper right',ncol=2)


fig.text(0.05, 0.5, r'$\rho(\varepsilon_{\rm KS}) $', va='center', rotation='vertical',fontsize=axis_fontsize)
plot(False,name='dos_'+phase+'.png')