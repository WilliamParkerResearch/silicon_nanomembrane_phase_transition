import numpy as np
#This is the data obtained from quantum esspresso
#pz-n-rrkjus_psl.0.1
#Diamond
total_energies_ecut_diamond = np.array([-90.93666425, -90.99836152, -91.00432167, -91.00489559, -91.00496606, -91.00504342, -91.00506127, -91.00508150, -91.00509117])
cutoff_energies_diamond = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90])
total_energies_kpoint_diamond = np.array([-90.99836152, -91.05527203, -91.06306773, -91.06447994, -91.06471799, -91.06478818, -91.06477768, -91.06481145, -91.06481099, -91.06480016, -91.06480297, -91.06479987, -91.06480366, -91.06480969, -91.06480398])
kpoints_diamond = np.array([2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16])
total_energies_strain_diamond = np.array([-90.90136971, -90.93771664, -90.96438048, -90.98312484, -90.99405063, -90.99836152, -90.99561333, -90.98728732, -90.97415241, -90.95572393, -90.93303237])
lattice_parameters_diamond = np.array([9.71938, 9.82169, 9.924, 10.0263, 10.1286, 10.2309, 10.3332, 10.4355, 10.5379, 10.6402, 10.7425])
celldm_1_diamond = 10.23092467
celldm_3_diamond = 1


#BetaSn
total_energies_ecut_BetaSn = np.array([-45.41870167, -45.45288626, -45.45662289, -45.45706090, -45.45708043, -45.45713454, -45.45712869, -45.45713772, -45.45714660])
cutoff_energies_BetaSn = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90])
total_energies_kpoint_BetaSn = np.array([-45.46381486, -45.37033540, -45.42209581, -45.47136946, -45.45788808, -45.44607494, -45.44904071, -45.45763130, -45.45270592, -45.45160494, -45.45316909, -45.45309711, -45.45357349, -45.45296173, -45.45257293])
kpoints_BetaSn = np.array([2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16])
total_energies_strain_all_BetaSn = np.array([-45.42380755, -45.43785525, -45.44784393, -45.45415028, -45.45712117, -45.45708043, -45.45432363, -45.44913395, -45.44175224, -45.43242864, -45.42137072])
lattice_parameters_all_BetaSn = np.array([9.00392, 9.0987, 9.19348, 9.28826, 9.38304, 9.47781, 9.57259, 9.66737, 9.76215, 9.85693, 9.9517])
total_energies_strain_shape_BetaSn = np.array([-45.45379728, -45.46608426, -45.47437200, -45.47902738, -45.48039400, -45.47879805, -45.47454245, -45.46788969, -45.45911409, -45.44844633, -45.43611428])
lattice_parameters_shape_BetaSn = np.array([9.00392, 9.0987, 9.19348, 9.28826, 9.38304, 9.47781, 9.57259, 9.66737, 9.76215, 9.85693, 9.9517])
celldm_1_BetaSn = 9.477779776
celldm_3_BetaSn = 0.48203609

# Unit conversions
total_energies_ecut_diamond = (13.601e-3/8)*total_energies_ecut_diamond                                     #(MeV/atom)
total_energies_kpoint_diamond = (13.601e-3/8)*total_energies_kpoint_diamond
total_energies_ecut_BetaSn = (13.601e-3/4)*total_energies_ecut_BetaSn
total_energies_kpoint_BetaSn = (13.601e-3/4)*total_energies_kpoint_BetaSn
total_energies_strain_diamond = (2.1798741e-18/8)*total_energies_strain_diamond                           #(Joules/atom)
total_energies_strain_all_BetaSn = (2.1798741e-18/4)*total_energies_strain_all_BetaSn
total_energies_strain_shape_BetaSn = (2.1798741e-18/4)*total_energies_strain_shape_BetaSn
cutoff_energies_diamond = 13.601e-3*cutoff_energies_diamond                                                     #(MeV)
cutoff_energies_BetaSn = 13.601e-3*cutoff_energies_BetaSn
lattice_parameters_diamond = 5.29177e-11*lattice_parameters_diamond                                           #(meters)
lattice_parameters_BetaSn = 5.29177e-11*lattice_parameters_all_BetaSn
celldm_1_diamond = 5.29177e-11*celldm_1_diamond
celldm_1_BetaSn = 5.29177e-11*celldm_1_BetaSn

