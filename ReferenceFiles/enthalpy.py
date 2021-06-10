import numpy as np


# H(V) = E + PV
#      = E(V) + P(V) V
def enthalpy_from_volume(parameters, volumes, eos='birch-murnaghan'):
    from equations_of_state import pressure_from_energy_equation_of_state
    if eos == 'birch-murnaghan':
        # H(V) = E(V) + (3 B0 / 2)((V0 / V)^(7/3) - (V0 / V)^(5/3)) * (1 + (3/4)(K0' - 4)((V0 / V)^(2/3) - 1))
        #       Sengupta et al., PRB 97, 235136 (2018) Eqns. 3 & 4
        from equations_of_state import birch_murnaghan
        pressures = pressure_from_energy_equation_of_state(parameters, volumes, eos=eos)
        energies = birch_murnaghan(parameters, volumes)
        enthalpies = energies + pressures * volumes
        return enthalpies
    elif eos == 'murnaghan':
        from equations_of_state import murnaghan
        pressures = pressure_from_energy_equation_of_state(parameters, volumes, eos=eos)
        energies = murnaghan(parameters, volumes)
        enthalpies = energies + pressures * volumes
        return enthalpies
    elif eos == 'vinet':
        from equations_of_state import vinet
        pressures = pressure_from_energy_equation_of_state(parameters, volumes, eos=eos)
        energies = vinet(parameters, volumes)
        enthalpies = energies + pressures * volumes
        return enthalpies
    else:
        print(f'No support yet for {eos} H(V)')
        return


def enthalpy_from_pressure(parameters, pressures, eos='murnaghan'):
    if eos == 'murnaghan':
        # H(P) = E0 + (K0 V0 /(K0' - 1)) * ( (1+(K0'/K0)*P)^(K0'-1)/K0' - 1)
        #           Shahi et al., PRB 97, 094111 (2018) Eqn. 3
        k_prime_minus_one = parameters[2] - 1.
        enthalpies = parameters[0] + \
                     (parameters[1] * parameters[3] / k_prime_minus_one) * \
                     (np.power(1. + parameters[2] / parameters[1] * pressures, k_prime_minus_one) / k_prime_minus_one)
        return enthalpies
    else:
        print(f'No support yet for {eos} H(P)')
        return
