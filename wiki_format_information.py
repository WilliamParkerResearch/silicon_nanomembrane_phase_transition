import numpy as np
from eos_information import *

# parameters = (E0, K0, K0', V0)
print('enertgy fit by munhagen')
print('')
print('|-')
print('| method')
print('| pseudoname')
print('|',round(1e30*fit_parameters_diamond[3], 2))
print('|',round(1e-9*fit_parameters_diamond[1], 1))
print('|',round(fit_parameters_diamond[2], 2))
print('|',round(1e30*tvol_diamond, 2))
print('|',round(1e30*fit_parameters_shape_BetaSn[3], 2))
print('|',round(celldm_3_BetaSn, 3))
print('|',round(1e-9*fit_parameters_shape_BetaSn[1], 1))
print('|',round(fit_parameters_shape_BetaSn[2], 2))
print('|',round(1e30*tvol_beta,2))
print('|',round(1e-9*tpressure, 2))
