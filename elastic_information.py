import numpy as np
import matplotlib.pyplot as plt
from importlib import import_module



def elastic_properties(N_ML,exchange_correlation = 'PBE',folder='elastic'):
    import numpy as np
    from importlib import import_module


    directoryofdata = 'DataFolder' + '.' + exchange_correlation + '.' + folder + '.' + 'Data_' + str(N_ML) + 'L'
    data = import_module(directoryofdata)

    strain_c_diamond = data.strain_c_diamond
    strain_a_diamond = data.strain_a_diamond
    strain_c_betasn = data.strain_c_betasn
    strain_a_betasn = data.strain_a_betasn
    stress_diamond = data.stress_diamond
    stress_betasn = data.stress_betasn

    poisson_d, b_d = -np.polyfit(strain_c_diamond,strain_a_diamond,1)
    poisson_b, b_b = -np.polyfit(strain_c_betasn,strain_a_betasn,1)

    youngs_d, by_d = -np.polyfit(strain_c_diamond,stress_diamond,1)
    youngs_b, by_b = -np.polyfit(strain_c_betasn,stress_betasn,1)
    return poisson_d,poisson_b,youngs_d,youngs_b


# p0 = elastic_properties(0)
# p1 = elastic_properties(1)
# p2 = elastic_properties(2)
# print(p0)
# print('----------------------------------------------------------------')
# print(p1)
# print((p1[0]-p0[0])/p0[0],(p1[1]-p0[1])/p0[1],(p1[2]-p0[2])/p0[2],(p1[3]-p0[3])/p0[3])
# print('----------------------------------------------------------------')
# print(p2)
# print((p2[0]-p0[0])/p0[0],(p2[1]-p0[1])/p0[1],(p2[2]-p0[2])/p0[2],(p2[3]-p0[3])/p0[3])
