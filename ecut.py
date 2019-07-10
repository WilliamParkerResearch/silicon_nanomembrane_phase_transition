import numpy as np
from scipy.interpolate import make_interp_spline
import matplotlib.pyplot as plt

etot = [-62.94878378, -62.95591968, -62.95659718, -62.95671392, -62.95673113, -62.95674196, -62.95675286, -62.95675978]
ecut = [20, 30, 40, 50, 60, 70, 80, 90]

xnew = np.linspace(ecut[0],ecut[-1],300)
spl = make_interp_spline(ecut, etot, k=3)
etot_smooth = spl(xnew)

plt.plot(xnew,etot_smooth)
plt.scatter(ecut, etot)
plt.title("Cutoff Energy vs Total Energy")
plt.xlabel("Cutoff Energy (Ry)")
plt.ylabel("Total Energy (Ry)")
plt.axhline(y=-62.9567, color='r')
plt.show()