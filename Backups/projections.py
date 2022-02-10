import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParamsDefault


data = np.loadtxt('DataFolder/PBE/bands/proj2mlp/Si.Fd-3m_2ML.k.pdos_tot')

k = np.unique(data[:,0])
e = np.unique(data[:,1])

dos = np.zeros([len(k), len(e)])

for i in range(len(data)):
    e_index = int(i % len(e))
    k_index = int(data[i][0] - 1)
    dos[k_index, e_index] = data[i][3]

plt.pcolormesh(k, e, dos.T)
plt.xticks([])
plt.ylabel('Energy (eV)')
plt.show()