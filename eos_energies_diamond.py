import matplotlib.pyplot as plt
from eos_information import *


# Use TeX fonts
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.sans-serif'] = "cmr10"


#   Diamond structure plot
fig7 = plt.figure()
plt.plot(volumes_diamond, fit_total_energies_strain_diamond*6.242e18 / n_atom_diamond)
plt.scatter(volumes_sim_diamond, total_energies_strain_diamond*6.242e18 / n_atom_diamond)
plt.title(r'Energy Equation of State for ' + structure_names[0] + ' ' + chemical_formula)
plt.xlabel(r'Volume (m$^3$/atom)')
plt.ylabel(r'Total Energy (eV/atom)')


plt.show()