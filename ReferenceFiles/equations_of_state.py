import numpy as np


def quadratic_fit(array):
    """
    Returns a NumPy array of values evaluated to an quadratic polynomial fit
    :param array:                 NumPy array(2, N) :: x- and y-values to be fit
    :return:                      NumPy quadratic_coefficients(3) :: coefficients of the fit quadratic polynomial
    """
    x_values = array[0, :]
    y_values = array[1, :]

    quadratic_coefficients = np.polyfit(x_values, y_values, 2)

    return quadratic_coefficients


def fit_eos(volumes, energies, quadratic_coefficients, eos='vinet', number_of_points=50):
    """
    Returns a NumPy array of values evaluated to an equation of state fit
    :param volumes:                 NumPy array(N) :: volumes (x-values) to be fit
    :param energies:                NumPy array(N) :: energies (y-values) to be fit
    :param quadratic_coefficients:  list(3) :: coefficients of the quadratic polynomial already fit to the data
    :param eos:                     str :: equation of state name ('murnaghan', 'birch-murnaghan', 'vinet')
    :param number_of_points:        int, optional :: number of points to evaluate fit function on
    :return:                        NumPy array(number_of_points) :: equation of state fit evaluated on grid,
    starting at volumes[0] and ending at volume[-1]
    """
    from scipy.optimize import curve_fit

    # Dictionary holding lambda functions from current module.
    lambda_dictionary = {
        'vinet': globals()['vinet'],
        'murnaghan': globals()['murnaghan'],
        'birch-murnaghan': globals()['birch_murnaghan']
    }

    # Get extremes of data and calculate range

    minimum_volume = np.amin(volumes)
    maximum_volume = np.amax(volumes)

    #   axis of symmetry: x = -b / 2a
    quadratic_axis_of_symmetry = -quadratic_coefficients[1]/(2*quadratic_coefficients[0])
    #   minimum: y = -b^2 / (4 a)  + c
    quadratic_minimum = -quadratic_coefficients[1]**2/(4*quadratic_coefficients[0]) + quadratic_coefficients[2]
    #   bulk modulus: K_0 = 2 * a / V_0 for E(V) = a*V^2 + b*V + E0
    quadratic_bulk_modulus = 2. * quadratic_coefficients[0] / quadratic_axis_of_symmetry

    bulk_modulus_derivative = 3.7

    # Get realistic equation of state fit
    initial_parameters = [quadratic_minimum, quadratic_bulk_modulus,
                          bulk_modulus_derivative, quadratic_axis_of_symmetry]

    eos_parameters, eos_covariances = curve_fit(lambda_dictionary[eos.lower()], volumes, energies, p0=initial_parameters)
    fit_curve_volumes = np.linspace(minimum_volume, maximum_volume, num=number_of_points)
    eos_fit_curve = lambda_dictionary[eos.lower()](fit_curve_volumes,
                                                   eos_parameters[0], eos_parameters[1], eos_parameters[2], eos_parameters[3])
    return eos_fit_curve, eos_parameters


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
