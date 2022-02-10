import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

cella = 10.33196444
n_atom = 32
percentage = np.array([-0.015,-0.010,-0.005,0.000,0.005,0.010,0.015,0.020])
energies = np.array([-373.4287204990,-373.4342942739,-373.4363166485,-373.2778819485,-373.4370053705,-373.4364161866,-373.4330358245,-373.4235991218])
percentagetest = np.array([-0.015,-0.010,-0.005,0.000,0.005,0.010,0.015])
energiestest = np.array([-373.4288179685,-373.3623961050,-373.3592542253,-373.2778584019,-373.3546551558,-373.3481922776,-373.3355118975])
# percentage = np.array([-0.015,-0.010,-0.005,0.005,0.010,0.015,0.020])
# energies = np.array([-373.4287204990,-373.4342942739,-373.4363166485,-373.4370053705,-373.4364161866,-373.4330358245,-373.4235991218])
# percentagetest = np.array([-0.015,-0.010,-0.005,0.005,0.010,0.015])
# energiestest = np.array([-373.4288179685,-373.3623961050,-373.3592542253,-373.3546551558,-373.3481922776,-373.3355118975])

p_max = np.amax(percentage)
p_min = np.amin(percentage)
n_val = 1000
ps = np.linspace(p_min, p_max, n_val)

energies_func = interpolate.interp1d(percentage, energies, kind='cubic')
energies_func2 = interpolate.interp1d(percentage, energies, kind='quadratic')

energies_fit = energies_func(ps)
energies_fit2 = energies_func2(ps)
min_per = ps[np.argwhere(energies_fit == np.amin(energies_fit))].flatten()
min_per2 = ps[np.argwhere(energies_fit2 == np.amin(energies_fit2))].flatten()

print(cella)
print(cella*(1-0.005))
print(cella*(1+min_per))
print(cella*(1+min_per2))
print(min_per)
print(min_per2)

plt.plot(ps,energies_fit2)
plt.scatter(percentage,energies)
plt.show()
plt.plot(percentage, energies)
plt.plot(percentagetest, energiestest)
plt.show()