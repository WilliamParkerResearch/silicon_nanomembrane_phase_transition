import numpy as np
import matplotlib.pyplot as plt

exchange_correlation = 'PBE'
N_ML = '3'

directoryofdata='DataFolder'+'.'+exchange_correlation+'.'+'bands'+'.'+'Data_Bands_'+N_ML+'L'
exec(f'from {directoryofdata} import *')

for i in range(nbands):
    exec(f'plt.plot(kpoints_diamond,bands_diamond_{i+1})')

plt.show()