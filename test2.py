import numpy as np
from eos_information import eos_properties
import matplotlib.pyplot as plt



# poisson_d, b_d = -np.polyfit(np.log(cellc_diamond), np.log(cella_diamond), 1)
# youngs_d, by_d = -np.polyfit(strain_c_diamond, stress_diamond, 1)
#
# youngs_b, by_b = -np.polyfit(strain_c_betasn, stress_betasn, 1)



# print(youngs_d*1e-9)

# plt.plot(strain_c_diamond,strain_a_diamond)
# plt.plot(strain_c_diamond,stress_diamond)
# plt.plot(strain_c_betasn,stress_betasn)
# plt.show()
# e00 = eos_properties(0)e0e = eos_properties(0,eos_folder='eosz')
e10 = eos_properties(1)

e1e = eos_properties(1,eos_folder='eosz')
e20 = eos_properties(2)
e2e = eos_properties(2,eos_folder='eosz')
# print(e0e)
# print(e00)
print(e1e)
print(e10)
print(e2e)
print(e20)

# eos_properties(0)[0])/eos_properties(0)[0]*100)