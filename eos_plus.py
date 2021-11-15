import scipy.optimize as sp
import numpy as np
import matplotlib.pyplot as plt

#           parameters = (E0, K0, K0', a0, c0)
nat = 32
# cella = (1e-10)*np.array([5.41295465409,5.43469342780,5.45104656750,5.46197050852,5.49477513620,5.54972288756,5.63296873087])
# cellc = (1e-10)*np.array([6.79196093800,6.76436063900,6.74292895667,6.73623870100,6.70329166833,6.64894243700,6.56114439700])
# energies = ((2.1798741e-18)/nat)*np.array([-1494.7768155813,-1494.7786394558,-1494.7789342564,-1494.7787075524,-1494.7757916786,-1494.7635292785,-1494.7285545278])

cella = (1e-10)*np.array([5.38569974917,5.44010075673,5.41295465409,5.43469342780,5.45104656750,5.46197050852,5.49477513620,5.54972288756,5.63296873087])
cellc = (1e-10)*np.array([6.87503430267,6.76216184833,6.79196093800,6.76436063900,6.74292895667,6.73623870100,6.70329166833,6.64894243700,6.56114439700])
energies = (2.1798741e-18/nat)*np.array([-1494.7727907596,-1494.7788520855,-1494.7768155813,-1494.7786394558,-1494.7789342564,-1494.7787075524,-1494.7757916786,-1494.7635292785,-1494.7285545278])
volumes = np.power(cella, 2) * cellc


def square_differences(p, a, b, y, f):
    return np.power(f(p, a, b) - y, 2).sum()


def mid(x):
    return int((len(x)-1)/2)


def eos(parameters, cella, cellc):
    volumes = np.power(cella,2)*cellc/8
    vol0 = (np.power(parameters[3],2)*parameters[4])/8
    k0pm1 = parameters[2] - 1.0  # K_0' - 1
    return parameters[0] + (parameters[1] * vol0 *
                            (((1.0 / (parameters[2] * k0pm1)) * np.power((volumes / vol0), (-k0pm1))) +
                             (volumes / (parameters[2] * vol0)) - (1.0 / k0pm1)))


initial_parameters = (energies[mid(energies)], 1e11, 3.5,
                              cella[mid(cella)],cellc[mid(cellc)])

fit_parameters = sp.fmin(square_differences, initial_parameters,
                                 args=(cella, cellc, energies, eos), maxiter=100000)

cellap = np.linspace(cella[0],cella[-1],1000)
cellcp = np.linspace(cellc[0],cellc[-1],1000)
energiesp = eos(fit_parameters,cellap,cellcp)

min_idx = np.argwhere(energiesp == np.amin(energiesp)).flatten()[0]
plt.plot(cella,cellc)
plt.scatter(cella,cellc)
plt.show()
# print(cellap[min_idx]*1.889e10)
# plt.plot(cellap*1.889,eos(fit_parameters,cellap,cellcp))
# plt.scatter(cella*1.889,energies)
# plt.show()
#
# print(fit_parameters[3]*1.889e10)