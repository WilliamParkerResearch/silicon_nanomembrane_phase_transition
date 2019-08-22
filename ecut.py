import numpy as np
from scipy.interpolate import make_interp_spline
import matplotlib.pyplot as plt

convergence_threshold = 7.35e-5     #this is = 1meV
etot = 13.601e-3/8*np.array([-62.94878378, -62.95591968, -62.95659718, -62.95671392, -62.95673113, -62.95674196, -62.95675286, -62.95675978])
ecut = 13.601e-3*np.array([20, 30, 40, 50, 60, 70, 80, 90])

plt.plot(ecut,etot, label='Total energy')
plt.scatter(ecut, etot)
plt.title("Cutoff Energy vs Total Energy for Diamond Structure of Si")
plt.xlabel("Cutoff Energy (meV)")
plt.ylabel("Total Energy (meV/atom)")
plt.axhline(etot[-1], color='r', label='convergence energy')
plt.axhline(etot[-1]-convergence_threshold, color='r', linestyle='--', label='convergence threshold')
plt.axhline(etot[-1]+convergence_threshold, color='r', linestyle='--')
plt.legend()
# plt.axis([0, ecut[-1],etot[-1]-10*convergence_threshold, etot[-1]+10*convergence_threshold])
plt.show()