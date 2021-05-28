import numpy as np
import matplotlib.pyplot as plt
from ReferenceFiles.unit_conversions import *

#This is the data obtained from quantum esspresso
#pbe-nl-kjpaw_psl.1.0.0

# SCF cycle convergence threshold
convergence_threshold = 1     # in Ry = 1 meV

# Simulation details
chemical_formula = 'Si'
structure_names = ['$Fd\overline{3}m$', '$I4_{1}/amd$']  # Hermann-Mauguin notation for diamond & beta-Sn structures
n_atom = 56


volumes_sim_diamond = cubic_meters_per_cubic_angstrom*np.array([24.68353997231,24.94336670804,25.20319344447,25.46302018090,25.72284691803,25.98267365447,26.24250039090,26.50232712804,26.76215386447,27.02198060090,27.28180733663])
total_energies_strain_diamond = (joules_per_rydberg/n_atom)*np.array([-2616.12698726,-2616.16247912,-2616.18980760,-2616.20916931,-2616.22069612,-2616.22458745,-2616.22072090,-2616.20955564,-2616.19123164,-2616.16591506,-2616.13380109])
volumes_sim_betasn = cubic_meters_per_cubic_angstrom*np.array([18.76740781056,18.96495947199,19.16251113341,19.36006279430,19.55761445518,19.75516611661,19.95271777803,20.15026943892,20.34782109980,20.54537276123,20.74292442265])
total_energies_strain_betasn = (joules_per_rydberg/n_atom)*np.array([-2615.26301732,-2615.29594153,-2615.32008931,-2615.33646593,-2615.34584412,-2615.34884839,-2615.34597449,-2615.33771504,-2615.32456611,-2615.30727438,-2615.28658154])


# plt.plot(volumes_sim_betasn,total_energies_strain_betasn)
# plt.plot(volumes_sim_diamond,total_energies_strain_diamond)
# plt.show()