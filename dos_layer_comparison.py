import numpy as np
import matplotlib.pyplot as plt
from ReferenceFiles.plot_formating import *
from scipy.signal import find_peaks
from importlib import import_module
from matplotlib.lines import Line2D

#Use TeX fonts
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.sans-serif'] = "cmr10"
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['axes.formatter.useoffset'] = False

adj_l = 1.2

exchange_correlation = 'PBE'
phase = 'diamond'
directoryofdata0 = 'DataFolder' + '.' + exchange_correlation + '.' + 'bands' + '.' + phase + '.' + 'Data_Bands_0L'
data0 = import_module(directoryofdata0)
nbands0 = data0.nbands
fermi_energy = data0.fermi_energy

bands0_energies = np.array([])
for i in np.arange(nbands0):
    exec(f'band_{i + 1} = data0.band_{i + 1}')
    exec(f'bands0_energies=np.append(bands0_energies,band_{i + 1})')

bands0_energies = np.sort(bands0_energies - fermi_energy)
band_gap_bound1_idx = np.argsort(np.abs(bands0_energies))[0]


if bands0_energies[band_gap_bound1_idx] <= 0:
    band_gap_bound2_idx = band_gap_bound1_idx+1
else:
    band_gap_bound2_idx = band_gap_bound1_idx-1

band_gap_bound1 = bands0_energies[band_gap_bound1_idx]
band_gap_bound2 = bands0_energies[band_gap_bound2_idx]

print('================================')
print(band_gap_bound2,band_gap_bound1)
print(band_gap_bound2-band_gap_bound1)
print('================================')


def dos_fermi(N_ML,phase = 'diamond',exchange_correlation = 'PBE'):
    from importlib import import_module
    N_ML = str(N_ML)
    print(N_ML)
    directoryofdata='DataFolder'+'.'+exchange_correlation+'.'+'bands'+'.'+phase+'.'+'Data_Bands_'+N_ML+'L'
    data = import_module(directoryofdata)

    dos_energies = data.dos_energies
    dos = data.dos
    fermi_energy = data.fermi_energy


    if N_ML == '0':
        nat = 8
    else:
        nat= float(N_ML)*8

    def integral(x, y):
        delta_x = np.diff(x)
        zero_idx = np.argsort(np.abs(x))[0]
        y_0 = y[zero_idx]
        x_0 = x[zero_idx]

        area_array = np.array([])
        for i in np.arange(len(delta_x)):
            idx_bounds = np.sort(np.array([i, zero_idx]))
            idxs = np.arange(idx_bounds[0], idx_bounds[-1] + 1, 1)
            if i < zero_idx:
                area = -np.sum(y[idxs] * delta_x[idxs])
            else:
                area = np.sum(y[idxs] * delta_x[idxs])
            area_array = np.append(area_array, area)

        return area_array



    integrated_dos = np.diff(dos_energies)*dos[1:]/nat
    integrated_dos2 = integral(dos_energies,dos/nat)


    minfe,maxfe = fermi_energy-2, fermi_energy+2
    frangeidxs = np.where(np.logical_and(dos_energies[:-1]>=minfe,dos_energies[:-1]<=maxfe))

    peaksr, _ = find_peaks(-integrated_dos[frangeidxs]+np.amax(integrated_dos[frangeidxs]), height=0)

    peaksrp = np.array([])
    for j in integrated_dos[frangeidxs][peaksr]:
        idx = np.where(integrated_dos == j)
        peaksrp = np.append(peaksrp,idx)
    peaksrp = peaksrp.astype(int)

    peaks, _ = find_peaks(-integrated_dos+np.amax(integrated_dos), height=0,width=10)

    if N_ML == '1':
        fermi_idx = peaks[np.argmin(np.abs(dos_energies[:-1][peaks]))]
    else:
        fermi_idx = peaks[np.argsort(integrated_dos[peaks])[0]]


    fermi_idxr = peaksrp[np.argsort(integrated_dos[peaksrp])[0]]

    fermi_energy_sim = dos_energies[:-1][fermi_idxr]
    p_0 = dos[np.where(dos_energies-fermi_energy_sim == 0)[0][0]]/nat
    p_0_og = dos[np.argsort(np.abs(dos_energies-fermi_energy))[0]]/nat

    # print('==============================')
    # print(fermi_energy)
    # print(fermi_energy_sim)


    ########################################################
    # area calculation
    bounds = np.where(np.logical_and(dos_energies - fermi_energy_sim >= band_gap_bound1,
                                     dos_energies - fermi_energy_sim <= band_gap_bound2))
    print(dos[bounds[0]]-fermi_energy)
    area_p = integrated_dos2[bounds[0][-1]] - integrated_dos2[bounds[0][0]]


    # # integrated plot
    # plt.figure(figsize=fs_m_12)
    # plt.title(str(N_ML))
    # plt.plot(dos_energies[:-1],integrated_dos)
    # plt.vlines(x=dos_energies[peaksr],ymin=0,ymax=integrated_dos[peaksr])
    # plt.vlines(x=dos_energies[fermi_idxr],ymin=0,ymax=integrated_dos[fermi_idxr],color='yellow')
    # plt.xlim(-10,10)
    # plt.show()

    # # # #Dos plot
    # plt.figure(figsize= fs_m_13)
    # plt.title(str(N_ML))
    # plt.plot(dos_energies - fermi_energy, dos/nat)
    # plt.fill(dos_energies - fermi_energy, dos/nat, color='goldenrod', alpha=0.5)
    # plt.vlines(x=dos_energies[fermi_idx]-fermi_energy,ymin=0,ymax=p_0)
    #
    # plt.vlines(x=dos_energies[bounds[0]]-fermi_energy,ymin=0,ymax=dos[bounds]/nat,color='black')
    #
    # plt.xlim(-5,5)




    return p_0_og


n_layers = np.array([0,1,2,3,4,5,6,7,8])

p_0_arr = np.array([])
for i in n_layers:
    p_0_arr = np.append(p_0_arr,dos_fermi(i))


plt.figure(figsize=fs_m_13)
plt.plot(n_layers[1:],p_0_arr[1:],color=rgbcode_diamond,linewidth=universal_linewidth,markersize=marker_size)
plt.scatter(n_layers[1:],p_0_arr[1:],marker=mark_d,color=rgbcode_diamond,linewidth=universal_linewidth,s=marker_size)

plt.hlines(y=p_0_arr[0],xmin=n_layers[1],xmax=n_layers[-1],color=adjust_lightness(rgbcode_diamond,adj_l),linestyle='dashed',linewidth=universal_linewidth)

plt.xlabel(r'$N_{\rm ML}$',fontsize=axis_fontsize)
plt.ylabel(r'$\rho$/atom',fontsize=axis_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.yticks(fontsize=tick_fontsize)

lines = [Line2D([0], [0], color=rgbcode_diamond_m,marker=mark_d,linewidth=universal_linewidth),
         Line2D([0], [0], color=adjust_lightness(rgbcode_diamond,adj_l),linestyle='dashed',linewidth=universal_linewidth)]
labels = [r'$N_{\rm ML}$','Bulk']

plt.legend(lines,labels,loc='upper right')


from layer_comparisons import cell_c_diamond


def second_axis_diamond():
    axes1 = plt.gca()
    axes2 = axes1.twiny()
    xticks2 = axes1.get_xticks()
    axes2.set_xlabel(r'$d ({\rm \AA})$',fontsize=axis_fontsize)
    axes2.tick_params(labelsize=tick_fontsize-2)

    xlabels2 = np.around(cell_c_diamond[(xticks2[1:-1]-1).astype(int)]*1e10,decimals=2)
    xlabels2 = np.append(xlabels2,'')
    xlabels2 = np.append('',xlabels2)

    x1_limi = axes1.get_xlim()[0]
    x1_limf = axes1.get_xlim()[1]
    nml_diff = xticks2[1:-1][-1]-xticks2[1:-1][0]
    d_diff = xlabels2[1:-1][-1].astype(float)-xlabels2[1:-1][0].astype(float)
    x2_limi = x1_limi*d_diff/nml_diff
    x2_limf = x1_limf*d_diff/nml_diff
    axes2.set_xlim(x2_limi,x2_limf)
    return

second_axis_diamond()

plot(True,'nml_fermi_dos_diamond')