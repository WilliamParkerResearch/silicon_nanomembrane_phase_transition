import numpy as np
from ReferenceFiles.unit_conversions import *
n_atom = 32
cella_diamond = (1e-10)*np.array([6.19052271733,6.18978853067,6.18704556533,6.18637493867,6.18637493867,6.18637493867,6.18637493867,6.18637493867,6.18565037467,6.18248366000,6.18145031067])
cellc_diamond = (1e-10)*np.array([5.17170288453,5.22614186240,5.28058084027,5.33501981813,5.38945879547,5.44389777333,5.49833675120,5.55277572853,5.60721470640,5.66165368427,5.71609266213])
stress_diamond = np.array([4841090677.774,3876079439.322,2907390566.578,1931052214.505,952654387.227,-25302123.935,-999286790.061,-1965327766.114,-2921071366.148,-3860486269.921,-4790633535.275])
strain_a_diamond = np.array([0.00067046998,0.00055179197,0.00010840382,0.00000000000,0.00000000000,0.00000000000,0.00000000000,0.00000000000,-0.00011712255,-0.00062900789,-0.00079604422])
strain_c_diamond = np.array([-0.05000000002,-0.04000000000,-0.02999999997,-0.01999999995,-0.01000000002,0.00000000000,0.01000000003,0.01999999995,0.02999999998,0.04000000000,0.05000000002])
total_energies_strain_diamond = (joules_per_rydberg/n_atom)*np.array([-1494.72712190,-1494.74605927,-1494.76059731,-1494.77090485,-1494.77702169,-1494.77900039,-1494.77691928,-1494.77094116,-1494.76113125,-1494.74767171,-1494.73060656])
volumes_sim_diamond = cella_diamond*cella_diamond*cellc_diamond/8
n_atom = 32
cella_betasn = (1e-10)*np.array([7.78905444400,7.76770619600,7.75614916600,7.74973021800,7.74982249200,7.75141470400,7.75434481800,7.75824410000,7.77635596800,7.78423078800,7.79009266600])
cellc_betasn = (1e-10)*np.array([4.94495097729,4.99700309290,5.04905520852,5.10110732387,5.15315943948,5.20521155510,5.25726367071,5.30931578632,5.36136790168,5.41342001729,5.46547213290])
stress_betasn = np.array([5508507749.245,4177204135.223,2963731923.945,1898100611.241,897048556.719,-18829487.580,-883956178.636,-1698037305.708,-2403113352.339,-3030517762.704,-3568923423.181])
strain_a_betasn = np.array([0.00485585425,0.00210174434,0.00061078683,-0.00021731336,-0.00020540921,0.00000000000,0.00037801022,0.00088105156,0.00321764026,0.00423356061,0.00498979392])
strain_c_betasn = np.array([-0.05000000001,-0.04000000000,-0.02999999999,-0.02000000002,-0.01000000001,0.00000000000,0.01000000001,0.02000000002,0.02999999999,0.04000000000,0.05000000001])
total_energies_strain_betasn = (joules_per_rydberg/n_atom)*np.array([-1494.36684143,-1494.38435797,-1494.39698514,-1494.40545457,-1494.41032818,-1494.41189467,-1494.41037109,-1494.40588925,-1494.39870156,-1494.38914426,-1494.37759102])
volumes_sim_betasn = cella_betasn*cella_betasn*cellc_betasn/8