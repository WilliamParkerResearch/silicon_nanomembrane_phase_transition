import matplotlib.pyplot as plt
from  eos_information import *

# Use TeX fonts
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.sans-serif'] = "cmr10"


#   Beta-Sin structure plot
fig8 = plt.figure()
plt.plot(volumes_BetaSn, fit_total_energies_strain_all_BetaSn*6.242e18 / n_atom_betasn)
plt.scatter(volumes_sim_BetaSn, total_energies_strain_all_BetaSn*6.242e18 / n_atom_betasn)
plt.plot(volumes_BetaSn, fit_total_energies_strain_shape_BetaSn * 6.242e18 / n_atom_betasn)
plt.scatter(volumes_sim_BetaSn, total_energies_strain_shape_BetaSn * 6.242e18 / n_atom_betasn)
plt.title(r'Energy Equation of State for ' + structure_names[1] + ' ' + chemical_formula)
plt.ylabel(r'Total Energy (eV/atom)')
plt.xlabel(r'Volume (m$^3$/atom)')


plt.show()