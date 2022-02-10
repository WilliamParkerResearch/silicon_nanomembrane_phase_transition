import matplotlib.pyplot as plt
from matplotlib import rcParamsDefault
import numpy as np
from ReferenceFiles.dos_functions import *
from ReferenceFiles.plot_formating import *
from ReferenceFiles.FunctionDefinitions import zpos_sort_idx, find_nearest,spline_prep
from ReferenceFiles.vasp_python_converter import vasp_data_modifier

exchange_correlation = 'PBE'
phase = 'diamond'
N_ML = '1'
order_idx = zpos_sort_idx(N_ML,phase)
zpos = vasp_data_modifier(N_ML,phase)[0][:,2][order_idx]-vasp_data_modifier(N_ML,phase)[0][:,2][order_idx][0]
nat = len(zpos)

if phase == 'diamond':
    prefix = 'Si.Fd-3m_'
    rgbcode_electric = rgbcode_diamond
    mark_e = mark_d
elif phase == 'betasn':
    prefix = 'Si.I4_1amd_'
    rgbcode_electric = rgbcode_betasn
    mark_e = mark_b

directoryofdata='DataFolder'+'.'+exchange_correlation+'.'+'bands'+'.'+phase+'.'+'Data_Bands_'+str(N_ML)+'L'
exec(f'from {directoryofdata} import fermi_energy')

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
    # t = [zpos[1]-(zpos_smod[1]-zpos_smod[0])/50 ,zpos[1]+(zpos_smod[3]-zpos_smod[1])/2, zpos[3]+(zpos_smod[4]-zpos_smod[3])/2]
    t = [(zpos_smod[0]+zpos_smod[-1])/2]
elif N_ML == '2':
    # t = [zpos_smod[1]+(zpos_smod[2]-zpos_smod[1])/4,zpos_smod[2],zpos_smod[4], zpos_smod[5], zpos_smod[-2]]
    t = [zpos_smod[2],zpos_smod[-2]]
elif N_ML == '3':
    # t = [zpos[1]+(zpos_smod[2]-zpos_smod[1])/4,zpos_smod[2],zpos_smod[5],zpos_smod[-2]]
    t = [zpos_smod[2],zpos_smod[-2]]
elif N_ML == '8':
    # t = [zpos[1]+(zpos_smod[2]-zpos_smod[1])/30,zpos_smod[3],zpos_smod[4]+(zpos_smod[5]-zpos_smod[4])/4,
         # zpos_smod[5]+(zpos_smod[6]-zpos_smod[5])/4,zpos_smod[6],zpos_smod[16],zpos_smod[-2]]
    t = [zpos_smod[5],zpos_smod[-2]]


dos_0_per_arr_sim = LSQUnivariateSpline(zpos_smod,dos_0_sarr,t,k=2)

fig = plt.figure(figsize=fs_s_12)

plt.plot(zpos_sim, dos_0_per_arr_sim(zpos_sim),color=rgbcode_electric,linewidth=universal_linewidth,zorder=1)
plt.scatter(zpos, dos_0_per_arr,color=rgbcode_electric,linewidth=universal_linewidth,marker=mark_e,zorder=2)


plt.ylim(0,0.31)
plt.gca().set_ylim(bottom=0)
plt.xlabel(r'$\mathit{z} ({\rm \AA})$',fontsize=axis_fontsize)
plt.ylabel(r'${\rho}_{atom}/{\rho}_{total}  $',fontsize=axis_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.yticks(fontsize=tick_fontsize)
plot(True,name='pdos_atom_contributios_'+phase+'_'+N_ML+'ml.png')
