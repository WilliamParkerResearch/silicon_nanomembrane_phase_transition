import numpy as np
from scipy.interpolate import make_interp_spline
import matplotlib.pyplot as plt

convergence_threshold = 7.35e-5     #this is = 1meV

etot = 13.601e-3/8*np.array([-62.27529336, -62.94878378, -63.00180153, -63.00875321, -63.01004222, -63.01028430, -63.01035681, -63.01034173, -63.01033394, -63.01035116])
kpoint = np.array([1,2,3,4,5,6,7,8,9,10])

plt.plot(kpoint,etot, label='total energy')
plt.scatter(kpoint, etot)
plt.title("Kpoints vs Total Energy of Diamond Structure of Si")
plt.xlabel("Number of Kpoints")
plt.ylabel("Total Energy (meV/atom)")
plt.axhline(etot[-1], color='r', label='convergence energy')
plt.axhline(etot[-1]-convergence_threshold, color='r', linestyle='--', label='convergence threshold')
plt.axhline(etot[-1]+convergence_threshold, color='r', linestyle='--')
plt.axis([0, kpoint[-1],etot[-1]-10*convergence_threshold, etot[-1]+10*convergence_threshold])
plt.legend()
plt.show()
