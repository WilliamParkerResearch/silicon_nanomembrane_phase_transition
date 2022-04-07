import matplotlib.pyplot as plt
from matplotlib import rcParamsDefault
import numpy as np
from ReferenceFiles.dos_functions import *
from ReferenceFiles.plot_formating import *
from ReferenceFiles.FunctionDefinitions import zpos_sort_idx


exchange_correlation = 'PBE'
phase = 'diamond'
N_ML = '22'
order_idx = zpos_sort_idx(N_ML,phase)
if N_ML == 0:
    nat = 8
else:
    nat = int(float(N_ML)*8)

if phase == 'diamond':
    prefix = 'Si.Fd-3m_'
    rgbcode_electric = rgbcode_diamond
    mark_e = mark_d
elif phase == 'betasn':
    prefix = 'Si.I4_1amd_'
    rgbcode_electric = rgbcode_betasn
    mark_e = mark_b

directoryofdata='DataFolder'+'.'+exchange_correlation+'.'+'bands'+'.'+phase+'.'+'Data_Bands_'+N_ML+'L'
exec(f'from {directoryofdata} import *')

pref = 'DataFolder/'+exchange_correlation+'/pdos/'+phase + '/'+str(N_ML)+'ML/'+prefix+str(N_ML)+'ML.k.pdos_'

et, pdost = pdos(pref+'tot','p')
norm_pdost = pdost/nbands
integrated_dos = np.diff(et)*norm_pdost[1:]



plt.plot(et[1:],integrated_dos)
plt.vlines(x=nat,ymin=np.amin(integrated_dos),ymax=np.amax(integrated_dos))
plt.show()