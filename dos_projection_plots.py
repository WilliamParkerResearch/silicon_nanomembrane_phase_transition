import matplotlib.pyplot as plt
from matplotlib import rcParamsDefault
import numpy as np
from ReferenceFiles.dos_functions import *
from ReferenceFiles.plot_formating import *
from ReferenceFiles.FunctionDefinitions import zpos_sort_idx


exchange_correlation = 'PBE'
phase = 'diamond'
N_ML = '1'
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
exec(f'from {directoryofdata} import *')

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


# atoms_dos = atom_pdos((np.arange(1,9,1)))
nel = nat*4


fig = plt.figure(figsize=fs_s_13)
plt.fill_between(et - fermi_energy,pdost/nbands,color=rgbcode_electric)
plt.plot(et - fermi_energy,pdost/nbands,color=adjust_lightness(rgbcode_electric,0.75),linewidth=universal_linewidth)

# plt.fill_between(et - fermi_energy,atoms_dos,color=adjust_lightness(rgbcode_electric,0.75))
# plt.plot(et - fermi_energy,atoms_dos,color=adjust_lightness(rgbcode_electric,0.60),linewidth=universal_linewidth)

plt.hlines(y=0,xmin=np.amin(et)-fermi_energy,xmax=np.amax(et)-fermi_energy,color=rgbcode_black,linewidth=universal_linewidth)
plt.xlim(-5,5)
plt.gca().set_ylim(bottom=0)
plt.xlabel(r'$\varepsilon_{\rm KS}$ (eV)',fontsize=axis_fontsize)
plt.ylabel(r'$\rho(\varepsilon_{\rm KS}) $',fontsize=axis_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.yticks([])
plot(True,name='dos_'+phase+'_'+N_ML+'ml.png')