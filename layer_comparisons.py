import numpy as np
import matplotlib.pyplot as plt
from eos_information import eos_properties
from ReferenceFiles.cohesive_energy_values import *
from ReferenceFiles.unit_conversions import *

n_layers = np.arange(0,8,1)
transition_pressures = np.array([])
transition_volumes_diamond = np.array([])
transition_volumes_betasn = np.array([])
volumes0_diamond = np.array([])
volumes0_betasn = np.array([])
energy0_diamond = np.array([])
energy0_betasn = np.array([])
bulk_mudulus0_diamond = np.array([])
bulk_mudulus0_betasn = np.array([])

for i in n_layers:
    if i == 0:
        t_pres0, vol_d0, vol_b0, fitp_d0, fitp_b0 = eos_properties(i)
        vol0_d0 = fitp_d0[3]
        energy0_d0 = fitp_d0[0]
        bulk_mod0_d0 = fitp_d0[1]
        vol0_b0 = fitp_b0[3]
        energy0_b0 = fitp_b0[0]
        bulk_mod0_b0 = fitp_b0[1]
    else:
        t_pres, vol_d, vol_b, fitp_d, fitp_b = eos_properties(i)
        transition_pressures = np.append(transition_pressures,t_pres)
        transition_volumes_diamond = np.append(transition_volumes_diamond,vol_d)
        transition_volumes_betasn = np.append(transition_volumes_betasn,vol_b)
        volumes0_diamond = np.append(volumes0_diamond,fitp_d[3])
        volumes0_betasn = np.append(volumes0_betasn,fitp_b[3])
        energy0_diamond = np.append(energy0_diamond,fitp_d[0])
        energy0_betasn = np.append(energy0_betasn,fitp_b[0])
        bulk_mudulus0_diamond = np.append(bulk_mudulus0_diamond,fitp_d[1])
        bulk_mudulus0_betasn = np.append(bulk_mudulus0_betasn,fitp_b[1])


cohesive_energy0_diamond = (energy0_diamond/8)-atom_energy_pbe
cohesive_energy0_d0 = (energy0_d0/8)-atom_energy_pbe
cohesive_energy0_betasn = (energy0_betasn/8)-atom_energy_pbe
cohesive_energy0_b0 = (energy0_b0/8)-atom_energy_pbe

plt.plot(n_layers[1:],-transition_pressures)
plt.plot(n_layers[1],-transition_pressures[0])
plt.hlines(y=-t_pres0,xmin=n_layers[1],xmax=n_layers[-1])
plt.hlines(y=0,xmin=n_layers[1],xmax=n_layers[-1],color='black')
plt.title('Transition Pressures')
plt.show()

plt.plot(n_layers[1:],transition_volumes_diamond)
plt.title('Transition Volumes Diamond')
plt.hlines(y=vol_d0,xmin=n_layers[1],xmax=n_layers[-1])
plt.show()

plt.plot(n_layers[1:],cohesive_energy0_diamond*eV_per_joule)
plt.title('Cohesive Energies Diamond')
plt.hlines(y=cohesive_energy0_d0*eV_per_joule,xmin=n_layers[1],xmax=n_layers[-1])
plt.show()

plt.plot(n_layers[1:],energy0_diamond)
plt.title('Ground State Energies Diamond')
plt.hlines(y=energy0_d0,xmin=n_layers[1],xmax=n_layers[-1])
plt.show()
#
plt.plot(n_layers[1:],volumes0_diamond)
plt.title('Ground State Volumes Diamond')
plt.hlines(y=vol0_d0,xmin=n_layers[1],xmax=n_layers[-1])
plt.show()

plt.plot(n_layers[1:],bulk_mudulus0_diamond)
plt.title('Bulk Modulus Diamond')
plt.hlines(y=bulk_mod0_d0,xmin=n_layers[1],xmax=n_layers[-1])
plt.show()

plt.plot(n_layers[1:],transition_volumes_betasn)
plt.hlines(y=vol_b0,xmin=n_layers[1],xmax=n_layers[-1])
plt.title('Transition Volumes BetaSn')
plt.show()

plt.plot(n_layers[1:],cohesive_energy0_betasn*eV_per_joule)
plt.title('Cohesive Energies BetaSn')
plt.hlines(y=cohesive_energy0_b0*eV_per_joule,xmin=n_layers[1],xmax=n_layers[-1])
plt.show()

plt.plot(n_layers[1:],energy0_betasn)
plt.title('Ground State Energies BetaSn')
plt.hlines(y=energy0_b0,xmin=n_layers[1],xmax=n_layers[-1])
plt.show()

plt.plot(n_layers[1:],volumes0_betasn)
plt.title('Ground State Volumes BetaSn')
plt.hlines(y=vol0_b0,xmin=n_layers[1],xmax=n_layers[-1])
plt.show()

plt.plot(n_layers[1:],bulk_mudulus0_betasn)
plt.title('Bulk Modulus BetaSn')
plt.hlines(y=bulk_mod0_b0,xmin=n_layers[1],xmax=n_layers[-1])
plt.show()

