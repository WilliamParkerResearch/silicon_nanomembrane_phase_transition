import numpy as np
#This is the data obtained from quantum esspresso
#Diamond
total_energies_ecut_diamond = np.array([-63.32390347, -63.34237499, -63.34571048, -63.34695836, -63.34719764, -63.34722498, -63.34722633, -63.34722701])
cutoff_energies_diamond = np.array([20, 30, 40, 50, 60, 70, 80, 90])
total_energies_kpoint_diamond = np.array([-62.63537678, -63.32390347, -63.37788374, -63.38498566, -63.38649362, -63.38668110, -63.38681263, -63.38676083, -63.38676166, -63.38679184])
kpoints_diamond = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
total_energies_strain_diamond = np.array([-62.93711345, -62.96873166, -62.99192902, -63.00739471, -63.01581990, -63.01783637, -63.01403802, -63.00497914, -62.99117832, -62.97312017, -62.95125931])
lattice_parameters_diamond = np.array([9.82927, 9.93274, 10.0362, 10.1397, 10.2431, 10.3466, 10.4501, 10.5535, 10.657, 10.7605, 10.8639])
celldm_1_diamond = 10.204655196
celldm_3_diamond = 1.000000000

#BetaSn
total_energies_ecut_BetaSn = np.array([-31.47687994, -31.62020167, -31.62960739, -31.63165071, -31.63237984, -31.63252938, -31.63254635, -31.63254733, -31.63254769])
cutoff_energies_BetaSn = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90])
total_energies_kpoint_BetaSn = np.array([-31.74799968, -31.63694751, -31.61916982, -31.62280005, -31.62766239, -31.63771738, -31.63154187, -31.63337729, -31.63157407, -31.63153220, -31.63303738, -31.63336299, -31.63268213, -31.63249175, -31.63221704])
kpoints_BetaSn = np.array([2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16])
total_energies_strain_all_BetaSn = np.array([-31.59634269, -31.60800569, -31.61575300, -31.62004952, -31.62124987, -31.61973952, -31.61556865, -31.60927137, -31.60079581, -31.59056224, -31.57867979])
lattice_parameters_all_BetaSn = np.array([8.91476, 9.0086, 9.10244, 9.19627, 9.29011, 9.38395, 9.47779, 9.57163, 9.66547, 9.75931, 9.85315])
total_energies_strain_z_BetaSn = np.array([-31.61873986, -31.62331915, -31.62541705, -31.62554737, -31.62371917, -31.62035202, -31.61549345, -31.60949341, -31.60236457, -31.59424356, -31.58526485])
lattice_parameters_z_BetaSn = np.array([8.91476, 9.0086, 9.10244, 9.19627, 9.29011, 9.38395, 9.47779, 9.57163, 9.66547, 9.75931, 9.85315])
celldm_1_BetaSn = 9.383953750
celldm_3_BetaSn = 0.506873898


# Unit conversions
total_energies_ecut_diamond = (13.601e-3/8)*total_energies_ecut_diamond                                     #(MeV/atom)
total_energies_kpoint_diamond = (13.601e-3/8)*total_energies_kpoint_diamond
total_energies_ecut_BetaSn = (13.601e-3/4)*total_energies_ecut_BetaSn
total_energies_kpoint_BetaSn = (13.601e-3/4)*total_energies_kpoint_BetaSn
total_energies_strain_diamond = (2.1798741e-18/8)*total_energies_strain_diamond                           #(Joules/atom)
total_energies_strain_all_BetaSn = (2.1798741e-18/4)*total_energies_strain_all_BetaSn
total_energies_strain_z_BetaSn = (2.1798741e-18/4)*total_energies_strain_z_BetaSn
cutoff_energies_diamond = 13.601e-3*cutoff_energies_diamond                                                     #(MeV)
cutoff_energies_BetaSn = 13.601e-3*cutoff_energies_BetaSn
lattice_parameters_diamond = 5.29177e-11*lattice_parameters_diamond                                           #(meters)
lattice_parameters_BetaSn = 5.29177e-11*lattice_parameters_all_BetaSn
celldm_1_diamond = 5.29177e-11*celldm_1_diamond
celldm_1_BetaSn = 5.29177e-11*celldm_1_BetaSn