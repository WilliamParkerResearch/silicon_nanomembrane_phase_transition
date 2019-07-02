import numpy as np
import scipy.optimize as sp
import matplotlib.pyplot as plt

bohr_to_meter = 5.29177e-11

lattice_parameters = np.array([9.81198545604, 9.91526951347, 10.01855357090, 10.12183762834, 10.22512168577,
                               10.32840574320, 10.43168980063, 10.53497385806, 10.63825791550, 10.74154197293,
                               10.84482603036])

volumes_sim = np.power(lattice_parameters * bohr_to_meter, 3)

# print(volumes_sim)

energies = 2.1798741e-19 * np.array([-62.93095553, -62.96400665, -62.98854005, -63.00524939, -63.01483051, -63.01791998,
                                     -63.01511204, -63.00696708, -62.99400720, -62.97671970, -62.95556260])

# print(energies)

initial_parameters = (np.min(energies), 1e11, 3.5, np.min(volumes_sim))

# parameters = (E0, K0, K0', V0)
# E(V) = E0 + (K0 V / K0')((V0/V)^K0' / (K0' - 1) + 1) - (K0 V0/(K0' - 1))

murnaghan = lambda p, v: p[0] + (p[1] * v / p[2]) * (np.power(p[3]/v, p[2]) / (p[2] - 1) + 1) - (p[1] * p[3]/(p[2] - 1))

square_differences = lambda p, x, y, f: np.power(f(p, x) - y, 2).sum()

fit_parameters = sp.fmin(square_differences, initial_parameters, args=(volumes_sim, energies, murnaghan), maxiter=100000)

print(fit_parameters)


volumes = np.linspace(1.3e-28, 2.0e-28, num=100)
fit_energies = murnaghan(fit_parameters, volumes)

plt.plot(volumes, fit_energies)
plt.scatter(volumes_sim, energies)
plt.axis([1.3e-28, 2.0e-28, -1.375e-17, -1.370e-17])

plt.show()