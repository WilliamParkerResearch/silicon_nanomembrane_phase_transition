import matplotlib.pyplot as plt
from matplotlib import rcParamsDefault
import numpy as np
from ReferenceFiles.FunctionDefinitions import spline_prep
from scipy.signal import find_peaks
from ReferenceFiles.plot_formating import *


n_layers = np.array([10])

def p_x(N_ML,phase='diamond',exchange_correlation='PBE'):
    N_ML = str(N_ML)
    if N_ML == '0':
        nat = 8
    else:
        nat = int(float(N_ML)*8)

    if phase == 'diamond':
        prefix = 'Si.Fd-3m_'
    elif phase == 'betasn':
        prefix = 'Si.I4_1amd_'


    directoryofdata='DataFolder'+'.'+exchange_correlation+'.'+'bands'+'.'+phase+'.'+'Data_Bands_'+N_ML+'L'
    exec(f'from {directoryofdata} import *')

    pref = 'DataFolder/'+exchange_correlation+'/pdos/'+phase + '/'+str(N_ML)+'ML/'+prefix+str(N_ML)+'ML.k.pdos_'
    fname = pref+'tot'





    def data_loader(fname):
        fid = open(fname, "r")
        data = fid.readlines()
        fid.close()

        energy = []
        pdos = []

        for row in range(len(data)):
            data_row = data[row]

            if (data_row[0][0] != '#' and len(data_row.split()) != 0):
                data_row = data_row[:-1].split()
                # print(data_row)
                energy.append(float(data_row[1]))
                pdos.append(float(data_row[3]))

        energy = np.asarray(energy)
        pdos = np.asarray(pdos)

        return energy, pdos



    energy, pdos_tot = data_loader(fname)

    energy_sorted = np.argsort(energy)

    et, pdost = spline_prep(energy[energy_sorted],pdos_tot[energy_sorted])

    #########################################################################################
    integrated_dos = np.diff(et)*pdost[1:]/nat

    peaks, _ = find_peaks(-integrated_dos+np.amax(integrated_dos), height=0,width=50)

    if N_ML == '1':
        fermi_idx = peaks[np.argmin(np.abs(et[:-1][peaks]))]
    else:
        fermi_idx = peaks[np.argsort(integrated_dos[peaks])[0]]
    # print(et[:-1][fermi_idx])


    fermi_energy_sim = et[:-1][fermi_idx]
    p_0 = pdost[np.where(et-fermi_energy_sim == 0)[0][0]]/nat
    #########################################################################################

    # make plots
    plt.figure(figsize = (8, 4))
    plt.plot(et-fermi_energy_sim, pdost/(nat), color='k', label='total')
    # plt.yticks([])
    plt.xlabel('Energy (eV)')
    plt.ylabel('DOS')

    plt.gca().set_ylim(bottom=0)

    plt.fill_between(et-fermi_energy_sim, pdost/nat, facecolor='k', alpha=0.25)
    plt.xlim(-5,5)
    plt.show()




    # # integrated plot
    # plt.figure(figsize=(8, 4))
    # plt.plot(et[:-1],integrated_dos)
    # plt.vlines(x=et[peaks],ymin=0,ymax=integrated_dos[peaks])
    # plt.vlines(x=et[fermi_idx],ymin=0,ymax=integrated_dos[fermi_idx],color='yellow')
    # plt.gca().set_ylim(bottom=0)
    #
    # plt.show()
    return p_0

p_0_arr = np.array([])
for i in n_layers:
    print(i)
    p_0_arr = np.append(p_0_arr,p_x(i))


plt.figure(figsize=fs_s_13)
plt.plot(n_layers[1:],p_0_arr[1:],color=rgbcode_diamond,linewidth=universal_linewidth,markersize=marker_size)
plt.scatter(n_layers[1:],p_0_arr[1:],marker=mark_d,color=rgbcode_diamond,linewidth=universal_linewidth,s=marker_size)

plt.hlines(y=p_0_arr[0],xmin=n_layers[1],xmax=n_layers[-1],color=rgbcode_diamond,linestyle='dashed',linewidth=universal_linewidth)

plt.xlabel('N-ML',fontsize=axis_fontsize)
plt.ylabel(r'$\rho(0) $',fontsize=axis_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.yticks(fontsize=tick_fontsize)
plot(True,'nml_fermi_dos_diamond')