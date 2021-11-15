import numpy as np
from ReferenceFiles.process_quantum_espresso_outputs import get_celldm, get_positions
from ReferenceFiles.pair_distribution import *

exchange_correlation = 'PBE'
phase = 'betasn'

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
h_n_ml = 7
# np.arange(l_n_ml, h_n_ml+1)
for i in np.arange(l_n_ml, h_n_ml+1):
    output_file = 'DataFolder'+'/'+exchange_correlation+'/'+'output_files'+'/'+prefix+str(i)+'ML.scf.out'
    lattice_parameter = get_celldm(output_file, verbose=verbose)
    lattice_parameter *= constants.value('Bohr radius')/constants.angstrom

    if maximum_distance <= 0.05:
        maximum_distance = np.ceil(lattice_parameter)

    if verbose:
        print('Maximum distance to plot g(r) = {} angstroms'.format(maximum_distance))
    if i == 0:
        atomic_positions = lattice_parameter * get_positions(output_file, verbose=verbose)

        rgbacode = rgbcode + (1,)
        plot_pdf(atomic_positions, density, minimum_distance=minimum_distance,maximum_distance=maximum_distance,verbose=verbose,linestyle='solid',color=rgbacode,linewidth=2.5)
    else:
        atomic_positions = lattice_parameter * multislab_builder(get_positions(output_file, verbose=verbose),3,3)
        minalpha = 0.4
        topalpha = 0.2
        pdiv = (i/(h_n_ml))
        alpha = minalpha*(1-pdiv)+pdiv*(1-topalpha)
        rgbacode = rgbcode + (alpha,)
        plot_pdf(atomic_positions, density, minimum_distance=minimum_distance,maximum_distance=maximum_distance,verbose=verbose,color=rgbacode,linewidth=1.25, dashes=(7,3))

if save_plot:
    plot_pdf_formating(atomic_positions, density,minimum_distance=minimum_distance, maximum_distance=maximum_distance, show_plot=False, verbose=verbose)
else:
    plot_pdf_formating(atomic_positions, density,minimum_distance=minimum_distance, maximum_distance=maximum_distance, show_plot=True, verbose=verbose)