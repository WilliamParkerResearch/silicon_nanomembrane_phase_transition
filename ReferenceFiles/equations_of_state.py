import numpy as np


def murnaghan(parameters, volumes):
    """
    Murnaghan equation of state: E(V) = E_0 + K_0 V_0 [ (1 / (K_0' (K_0' - 1))) (V / V_0)^(-(K_0' - 1)) +
                                                        (1 / K_0') (V / V_0) -

                                                        (1 / (K_0' - 1)) ]

    :param parameters: list of equation of state parameters  -  equilibrium energy (E_0),
                                                                bulk modulus (K_0),
                                                                bulk modulus pressure derivative (K_0'),
                                                                equilibrium volume (V_0)
    :param volumes: NumPy array of volumes per atom
    :return: NumPy array of Murnaghan equation of state values at input volumes
    """
    k0pm1 = parameters[2] - 1.0  # K_0' - 1
    return parameters[0] + (parameters[1] * parameters[3] *
                            (((1.0 / (parameters[2] * k0pm1)) * np.power((volumes / parameters[3]), (-k0pm1))) +
                             (volumes / (parameters[2] * parameters[3])) - (1.0 / k0pm1)))


def birch_murnaghan(parameters, volumes):
    """
    Birch-Murnaghan equation of state: E(V) = E_0 + (9/16) K_0 V_0 {[ (V / V_0)^(-(2/3)) - 1 ]^3 K_0' +
                                                                    [ (V / V_0)^(-(2/3)) - 1]^2 *
                                                                    [ 6 - 4 (V / V_0)^(-(2/3)) ]}
    :param parameters: list of equation of state parameters     equilibrium energy (E_0),
                                                                bulk modulus (K_0),
                                                                bulk modulus pressure derivative (K_0'),
                                                                equilibrium volume (V_0)
    :param volumes: NumPy array of volumes per atom
    :return: NumPy array of the Birch-Murnaghan equation of state values at input volumes
    """
    reduced_volume_area = np.power(volumes / parameters[3], -2. / 3.)
    return parameters[0] + (9. * parameters[1] * parameters[3] / 16.) * (
            np.power(reduced_volume_area - 1., 3.) * parameters[2] +
            np.power(reduced_volume_area - 1., 2.) * (6. - 4. * reduced_volume_area))


def vinet(parameters, volumes):
    """
    Vinet equation of state: E(V) = E_0 + (2 K_0 V_0 / (K_0' - 1)^2) *
                                        - {2 - [5 + 3 (V / V_0)^(1/3) (K_0' - 1) - 3 K_0']
                                               exp(- (3/2) (K_0' - 1) (1 - (V / V_0)^(1/3))})
    :param parameters: list of equation of state parameters     equilibrium energy (E_0),
                                                                bulk modulus (K_0),
                                                                bulk modulus pressure derivative (K_0'),
                                                                equilibrium volume (V_0)
    :param volumes: NumPy array of volumes per atom
    :return: NumPy array of the Vinet equation of state values at input volumes
    """
    k0pm1 = parameters[2] - 1  # K_0' - 1
    eta = 1.5 * k0pm1  # (3/2) (K_0' - 1)
    reduced_volume_lengths = np.power(volumes / parameters[3], 1. / 3.)  # x = (V/V_0)^(1/3)

    # Hama and Suito 1996, Eqn 16
    vinet_eos = parameters[0] + \
                9. * parameters[1] * parameters[3] / eta**2 * \
                (1. - (1. - eta*(1. - reduced_volume_lengths)) * np.exp(eta * (1. - reduced_volume_lengths)))

    return vinet_eos


def pressure_from_energy_equation_of_state(parameters, volumes, eos='birch-murnaghan'):
    if eos == 'birch-murnaghan':
        reduced_volumes = parameters[3] / volumes  # V_0 / V = v
        # P(V) = (3/2) K_0 (v^(7/3) - v^(5/3)) (1 + (3/4)(K_0' - 4)(v^(2/3) -1)
        pressures = 3 * parameters[1] / 2 * \
                    (np.power(reduced_volumes, 7 / 3) - np.power(reduced_volumes, 5 / 3)) * \
                    (1. + (3 / 4) * (parameters[2] - 4.) * (np.power(reduced_volumes, 2 / 3) - 1.))
        return pressures
    elif eos == 'murnaghan':
        # P(V) = (K_0 / K_0') ((V_0/V)^K_0' - 1)
        pressures = (parameters[1] / parameters[2]) * (np.power(parameters[3]/volumes, parameters[2]) - 1.)
        return pressures
    elif eos == 'vinet':
        reduced_volume_lengths = np.power(volumes / parameters[3], 1./3.)  # v = V/V_0
        # P(V) = 3 B0 (1 - v) v^-2 e^((3/2)(B0' -1)(1-v))
        pressures = 3. * parameters[1] * (1. - reduced_volume_lengths) * reduced_volume_lengths**-2 * \
                    np.exp((3./2) * (parameters[2] - 1.) * (1. - reduced_volume_lengths))
        return pressures
    else:
        print('No P(V) implemented yet for {}'.format(eos))
        return
