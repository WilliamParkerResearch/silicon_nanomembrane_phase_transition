import numpy as np
from ReferenceFiles.unit_conversions import *
n_atom = 24
cella_diamond = (1e-10)*np.array([5.41203821054,5.41203821054,5.41203821054,5.41203821054,5.41203821054,5.41203821054,5.41203821054,5.41203821054,5.41203821054,5.41203821054,5.41203821054])
cellb_diamond = (1e-10)*np.array([5.41203821054,5.41203821054,5.41203821054,5.41203821054,5.41203821054,5.41203821054,5.41203821054,5.41203821054,5.41203821054,5.41203821054,5.41203821054])
cellc_diamond = (1e-10)*np.array([5.03450838545454545454,5.08750321090909090909,5.14049803563636363636,5.19349286109090909090,5.24648768581818181818,5.29948251127272727272,5.35247733672727272727,5.40547216145454545454,5.45846698690909090909,5.51146181163636363636,5.56445663709090909090])
total_energies_strain_diamond = (joules_per_rydberg/n_atom)*np.array([-188.40581464,-188.42044794,-188.43170039,-188.43945291,-188.44450221,-188.44602339,-188.44443488,-188.43994921,-188.43231382,-188.42200028,-188.40853625])
volumes_sim_diamond = cella_diamond*cellb_diamond*cellc_diamond/8
n_atom = 24
cella_betasn = (1e-10)*np.array([4.77033373765,4.77033373765,4.77033373765,4.77033373765,4.77033373765,4.77033373765,4.77033373765,4.77033373765,4.77033373765,4.77033373765])
cellb_betasn = (1e-10)*np.array([4.77033373765,4.77033373765,4.77033373765,4.77033373765,4.77033373765,4.77033373765,4.77033373765,4.77033373765,4.77033373765,4.77033373765])
cellc_betasn = (1e-10)*np.array([4.93686880382608695652,4.98883584417391304347,5.04080288417391304347,5.09276992417391304347,5.14473696417391304347,5.19670400417391304347,5.24867104417391304347,5.30063808417391304347,5.35260512417391304347,5.40457216417391304347])
total_energies_strain_betasn = (joules_per_rydberg/n_atom)*np.array([-187.96703507,-187.98082734,-187.99103007,-187.99735054,-188.00102449,-188.00222632,-188.00082500,-187.99775329,-187.99245643,-187.98552183])
volumes_sim_betasn = cella_betasn*cellb_betasn*cellc_betasn/8