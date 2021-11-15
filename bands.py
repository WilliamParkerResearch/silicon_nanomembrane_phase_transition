import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

exchange_correlation = 'PBE'
phase = 'diamond'
N_ML = '2'

directoryofdata='DataFolder'+'.'+exchange_correlation+'.'+'bands'+'.'+phase+'.'+'Data_Bands_'+N_ML+'L'
exec(f'from {directoryofdata} import *')


fig = plt.figure()
gs = fig.add_gridspec(nrows=1, ncols=5)
ax1 = fig.add_subplot(gs[0,:4])
for i in range(nbands):
    exec(f'plt.plot(kpoints,band_{i+1}-fermi_energy,color="goldenrod")')
plt.hlines(y=0,xmin=kpoints[0],xmax=kpoints[-1],color='black',linestyle='dotted')
plt.xticks(kpoint_idx,kpoint_names)
ax1.set_ylim(-5,5)
ax1.set_xlim(kpoints[0],kpoints[-1])
ax1.set_ylabel(r'$\varepsilon_{\rm KS}$ (eV)')
ax2 = fig.add_subplot(gs[0,4:])
plt.plot(dos,dos_energies-fermi_energy)
ax2.fill(dos, dos_energies - fermi_energy, color='goldenrod', alpha=0.5)
ax2.set_ylim(-5,5)
ax2.set_xlim(0,np.amax(dos))
ax2.yaxis.tick_right()
plt.show()