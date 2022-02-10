import numpy as np
import matplotlib.pyplot as plt
from eos_information import eos_properties

N=6
n_cell = np.arange(1,N,1)
transition_pressures = np.array([])
transition_volumes_diamond = np.array([])
transition_volumes_betasn = np.array([])

t_pres0N, vol_d0N, vol_b0N,par_d0N,par_b0N = eos_properties(N)

for i in n_cell:
    t_pres, vol_d, vol_b, par_d, par_b= eos_properties(N,OF=True,n_etype=i,etype='ce')
    transition_pressures = np.append(transition_pressures,t_pres)
    transition_volumes_diamond = np.append(transition_volumes_diamond,vol_d)
    transition_volumes_betasn = np.append(transition_volumes_betasn,vol_b)


plt.plot(n_cell,-transition_pressures)
plt.hlines(y=-t_pres0N,xmin=n_cell[0],xmax=n_cell[-1],linestyle='dashed')
plt.title('Transition Pressures')
plt.xlabel('cell in slab from bottom to top')
plt.show()

plt.plot(n_cell,transition_volumes_diamond)
plt.hlines(y=vol_d0N,xmin=n_cell[0],xmax=n_cell[-1],linestyle='dashed')
plt.title('Transition Volumes Diamond')
plt.xlabel('cell in slab from bottom to top')
plt.show()

plt.plot(n_cell,transition_volumes_betasn)
plt.hlines(y=vol_b0N,xmin=n_cell[0],xmax=n_cell[-1],linestyle='dashed')
plt.title('Transition Volumes BetaSn')
plt.xlabel('cell in slab from bottom to top')
plt.show()