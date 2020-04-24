import numpy as np

# This is the data obtained from quantum esspresso
# pz-n-rrkjus_psl.0.1

n_atom_diamond = 8
n_atom_betasn = 4

# Diamond
total_energies_ecut_diamond = np.array([-90.99873344, -91.06316669, -91.06946699, -91.07005680, -91.07012621, -91.07020444, -91.07022241, -91.07024270, -91.07025222])
cutoff_energies_diamond = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90])
total_energies_kpoint_diamond = np.array([-90.99775717, -91.05532974, -91.06316669, -91.06455147, -91.06480878, -91.06487732, -91.06489348, -91.06491916, -91.06490842, -91.06490603, -91.06489977, -91.06490971, -91.06490218, -91.06490645, -91.06490465])
kpoints_diamond = np.array([2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16])
total_energies_strain_diamond = np.array([-90.96738411, -91.00330350, -91.02988291, -91.04849901, -91.05929854, -91.06316669, -91.06078685, -91.05271284, -91.03941758, -91.02142347, -90.99925471])
lattice_parameters_diamond = np.array([9.69089, 9.7929, 9.89491, 9.99692, 10.0989, 10.2009, 10.3029, 10.405, 10.507, 10.609, 10.711])
celldm_1_diamond = 10.20094043
celldm_3_diamond = 1

# BetaSn
total_energies_ecut_BetaSn = np.array([-45.43565043, -45.47025292, -45.47413256, -45.47452449, -45.47456329, -45.47461178, -45.47461769, -45.47463284, -45.47463542])
cutoff_energies_BetaSn = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90])
total_energies_kpoint_BetaSn = np.array([-45.58596237, -45.48359686, -45.46202177, -45.46382037, -45.46750173, -45.47964456, -45.47330140, -45.47526333, -45.47343472, -45.47320680, -45.47470593, -45.47524001, -45.47454320, -45.47456329, -45.47396820])
kpoints_BetaSn = np.array([2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16])
total_energies_strain_all_BetaSn = np.array([-45.43131918, -45.44780124, -45.45999767, -45.46828445, -45.47302502, -45.47456329, -45.47319140, -45.46922246, -45.46293770, -45.45457158, -45.44436659])
lattice_parameters_all_BetaSn = np.array([8.51788, 8.60755, 8.69721, 8.78687, 8.87653, 8.96619, 9.05586, 9.14552, 9.23518, 9.32484, 9.4145])
total_energies_strain_shape_BetaSn = np.array([-45.43183349, -45.44811691, -45.46015898, -45.46834415, -45.47303218, -45.47456851, -45.47323927, -45.46934743, -45.46315979, -45.45492431, -45.44489543])
lattice_parameters_shape_BetaSn = np.array([8.51788, 8.60755, 8.69721, 8.78687, 8.87653, 8.96619, 9.05586, 9.14552, 9.23518, 9.32484, 9.4145])
celldm_1_BetaSn = 8.96619465
celldm_3_BetaSn = 0.554456931

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
