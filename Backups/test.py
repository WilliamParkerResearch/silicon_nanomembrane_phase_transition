import numpy as np
from ReferenceFiles.process_quantum_espresso_outputs import get_celldm, get_positions
from ReferenceFiles.pair_distribution import *

exchange_correlation = 'PBE'
phase = 'diamond'

if phase == 'diamond':
    prefix = 'Si.Fd-3m_'
    density = 2300
    rgbcode = (60/255, 80/255, 155/255)
if phase == 'betasn':
    prefix = 'Si.I4_1amd_'
    density = 3370
    rgbcode = (156/255, 61/255, 61/255)


verbose = False
save_plot = False
minimum_distance = 1.5
maximum_distance = 6
l_n_ml = 0
h_n_ml = 1
square_l =3

print('-----------------------------------------------------------------------------------')

output_file = 'DataFolder'+'/'+exchange_correlation+'/'+'output_files'+'/'+prefix+str(1)+'ML.scf.out'
lattice_parameter = get_celldm(output_file, verbose=verbose)
lattice_parameter *= constants.value('Bohr radius')/constants.angstrom
print(lattice_parameter)

atomic_positions = get_positions(output_file, verbose=verbose)
def multislab_builder(positions,lx=3,ly=3):
    one_xarray = np.ones((len(positions),len(positions[0])))
    one_xarray[0:,1:]=0
    tmp_xarray = 1*positions
    empty_array = np.empty((0,3))
    empty_array=np.append(empty_array,tmp_xarray,axis=0)

    for i in np.arange(1,lx+1):
        tmp_xarray=tmp_xarray+one_xarray
        empty_array=np.append(empty_array,tmp_xarray,axis=0)

    one_yarray = np.ones((len(empty_array),len(empty_array[0])))
    one_yarray[0:,(0,-1)]=0
    tmp_yarray = 1*empty_array

    for i in np.arange(1,ly+1):
        tmp_yarray=tmp_yarray+one_yarray
        empty_array=np.append(empty_array,tmp_yarray,axis=0)
    return empty_array
print(multislab_builder(atomic_positions))