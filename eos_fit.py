import numpy as np
import scipy.optimize as sp
import matplotlib.pyplot as plt


A =8
B = 1
bohr_to_meter = 5.29177e-11

lattice_parameters = np.array([9.81198545604, 9.91526951347, 10.01855357090, 10.12183762834, 10.22512168577,
                               10.32840574320, 10.43168980063, 10.53497385806, 10.63825791550, 10.74154197293,
                               10.84482603036])
volumes_sim = (1/A)*np.power(lattice_parameters * bohr_to_meter, 3)

# print(volumes_sim)

energies = (2.1798741e-18/A) * np.array([-62.93095553, -62.96400665, -62.98854005, -63.00524939, -63.01483051, -63.01791998,
                                     -63.01511204, -63.00696708, -62.99400720, -62.97671970, -62.95556260])

# print(energies)

initial_parameters = (np.min(energies), 1e11, 3.5, np.min(volumes_sim))

# parameters = (E0, K0, K0', V0)
# E(V) = E0 + (K0 V / K0')((V0/V)^K0' / (K0' - 1) + 1) - (K0 V0/(K0' - 1))

def murnaghan(p, v):
    kk = p[2]-1.0
    # return p[0] + (p[1] * v / p[2]) * (np.power(p[3]/v, p[2]) / (p[2] - 1) + 1) - (p[1] * p[3]/(p[2] - 1))
    return (p[0] + (p[1]*p[3]*(((1.0/(p[2]*kk))*np.power((v/p[3]), (-kk)))+(v/(p[2]*p[3]))-(1.0/kk))))


def square_differences(p, x, y, f):
    return np.power(f(p, x) - y, 2).sum()

fit_parameters = sp.fmin(square_differences, initial_parameters, args=(volumes_sim, energies, murnaghan), maxiter=100000)

volumes = np.linspace(1.3e-28, 2.0e-28, num=100)
fit_energies = murnaghan(fit_parameters, volumes)

plt.plot(volumes, fit_energies*6.242e18/8)
plt.scatter(volumes_sim, energies*6.242e18/8)
# plt.axis([1.3e-28, 2.0e-28, -1.375e-17*6.242e18, -1.370e-17*6.242e18])
plt.title('Equation of State of Diamond Structure')
plt.ylabel('Total Energy (eV/atom)')
plt.xlabel('Volume (m^3)')


# plt.show()
# parameters = (E0, K0, K0', V0)
fit_parameters_conv = [6.2415093433e18*fit_parameters[0]/B, 1e-9*fit_parameters[1], fit_parameters[2], fit_parameters[3]/B]

print(fit_parameters)