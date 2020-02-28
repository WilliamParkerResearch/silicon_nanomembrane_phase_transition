import matplotlib.pyplot as plt
from eos_information import *


# Use TeX fonts
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.sans-serif'] = "cmr10"


# Two-phase plot
fig9 = plt.figure()
plt.plot(volume/cubic_meters_per_cubic_angstrom, (murnaghan(fit_parameters_shape_BetaSn, tvol_beta)-(volume-tvol_beta)*tpressure)*6.242e18,color='y')
plt.plot(volumes_diamond/cubic_meters_per_cubic_angstrom, fit_total_energies_strain_diamond*6.242e18,color='slateblue')
plt.scatter(volumes_sim_diamond/cubic_meters_per_cubic_angstrom, total_energies_strain_diamond*6.242e18,label=r'Diamond',color='slateblue',marker='s')
plt.plot(volumes_BetaSn/cubic_meters_per_cubic_angstrom, fit_total_energies_strain_shape_BetaSn*6.242e18,color='firebrick')
plt.scatter(volumes_sim_BetaSn/cubic_meters_per_cubic_angstrom, total_energies_strain_shape_BetaSn*6.242e18, label=r'BetaSn',color='firebrick',marker='v')
plt.yscale(r'linear')
plt.title(r'Energy Equation of State for ' + structure_names[0] + ' and ' + structure_names[1] + ' ' +
               chemical_formula)
plt.ylabel(r'Total Energy (eV/atom)')
plt.xlabel(r'Volume (\r{A}$^3$/atom)')
plt.legend()


plt.show()