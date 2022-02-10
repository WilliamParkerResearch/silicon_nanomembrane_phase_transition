import numpy as np
import matplotlib.pyplot as plt
from importlib import import_module


exchange_correlation = 'PBE'
phase = 'diamond'
directoryofdata0 = 'DataFolder' + '.' + exchange_correlation + '.' + 'bands' + '.' + phase + '.' + 'Data_Bands_10L'
data = import_module(directoryofdata0)

dos = data.dos
dos_energies = data.dos_energies
fermi_energy = data.fermi_energy

plt.plot(dos_energies-fermi_energy,dos)
plt.xlim(-5,5)
plt.show()