#!/usr/bin/env python3
import numpy as np
from scipy import constants


def pair_distribution(distances, radial_distributions, mass_density):
    """
              g(r) = G(r) / (4 pi r^2 rho)               [units = kg^-1]
              G(r) = radial distribution function        [units = angstrom^-1]
              r    = interatomic distances               [units = angstrom]
              rho  = mass density of reference structure [units = kg m^-3]
    """
    distances = distances * constants.angstrom                        # convert distances in angstrom to meter
    radial_distributions = radial_distributions / constants.angstrom  # convert rdf in angstrom^-1 to m^-1
    pdf = radial_distributions / (4 * np.pi * distances**2 * mass_density)
    return pdf * constants.angstrom**3


def radial_distribution(distances, number_of_atoms,
                        delta_width=0.015, minimum_distance=0., maximum_distance=0., number_of_points=10000,
                        verbose=False):
    """
              G(r)     = (1/N) sum_ij gaussian(r - |R_i - R_j|) [units = angstrom^-1]
              N        = number of atoms
              R_i, R_j = position of atoms i and j              [units = angstrom]
    """
    gaussian_value_threshold = 1.e-4
    if verbose:
        print('Calculating radial distribution...')

    if minimum_distance == 0. and maximum_distance == 0.:
        distance_values = np.linspace(np.amin(distances), np.amax(distances), num=number_of_points)
    else:
        distance_values = np.linspace(minimum_distance, maximum_distance, num=number_of_points)
    if verbose:
        print('\tSeparation values array shape = {}'.format(distance_values.shape))
        print('\tSeparation range = [{}, {}]'.format(np.amin(distance_values), np.amax(distance_values)))

    rdf = np.zeros(len(distance_values))
    for distance in distances:
        if verbose:
            print('For interatomic distance {}:'.format(distance))
        gaussian_values = gaussian(distance_values, mean=distance, width=delta_width)
        if np.amax(gaussian_values) < gaussian_value_threshold:
            gaussian_maximum = 1.
        else:
            gaussian_maximum = np.amax(gaussian_values)
        if verbose:
            print('\tGaussian scale = {}'.format(gaussian_maximum))
            print('\tGaussian values dimensions = {}'.format(gaussian_values.shape))
        rdf = rdf + gaussian_values / gaussian_maximum

    rdf = rdf / number_of_atoms

    return rdf, distance_values


def gaussian(x, mean=0., width=1.):
    """gaussian(x) = e^(-((x-mu)/sigma)^2 / 2) / (sqrt(2 pi) sigma)"""
    return np.exp(-0.5 * (x - mean)**2 / width**2)  / (np.sqrt(2*np.pi) * width)


def interatomic_distances(positions):
    """atomic positions should be numpy array"""
    distances = []
    atomic_index = 1
    for atom_position in positions:
        for second_atom_position in positions[atomic_index:]:
            distances.append(calculate_distance(atom_position, second_atom_position))
        atomic_index += 1
    return np.array(distances)


def calculate_distance(position_one, position_two):
    """position_one and position_two should be numpy arrays"""
    separation = position_one - position_two
    distance = np.sqrt(np.sum(separation ** 2))
    return distance


def plot_pdf(positions, mass_density,minimum_distance=0.05, maximum_distance=0.05, verbose=False,zorder=1,color='grey',linestyle='dashed',linewidth=1,dashes=(3,2)):
    distances = interatomic_distances(positions)
    # distances = distances[(distances >= minimum_distance) & (distances <= maximum_distance)]

    if verbose:
        print('Interatomic distances range = [{}, {}]'.format(np.amin(distances), np.amax(distances)))
        print('Interatomic distances array shape = {}'.format(distances.shape))
        print('Number of atoms to calculate = {}'.format(len(positions)))
    radial_distribution_values, distance_values = radial_distribution(distances, len(positions),
                                                                      minimum_distance=minimum_distance,
                                                                      maximum_distance=maximum_distance,
                                                                      verbose=verbose)
    pair_distribution_values = pair_distribution(distance_values, radial_distribution_values, mass_density)
    import matplotlib.pyplot as plt
    if linestyle == 'dashed' or linestyle == '--':
        plt.plot(distance_values, pair_distribution_values,zorder=zorder,color=color,linestyle=linestyle,linewidth=linewidth,dashes=dashes)
    else:
        plt.plot(distance_values, pair_distribution_values,zorder=zorder,color=color,linestyle=linestyle,linewidth=linewidth)

    return


def plot_pdf_formating(positions, mass_density,minimum_distance=0.05, maximum_distance=0.05, show_plot=True, verbose=False, filename='pdf.png', label=''):
    distances = interatomic_distances(positions)
    # distances = distances[(distances >= minimum_distance) & (distances <= maximum_distance)]

    # if verbose:
    #     print('Interatomic distances range = [{}, {}]'.format(np.amin(distances), np.amax(distances)))
    #     print('Interatomic distances array shape = {}'.format(distances.shape))
    #     print('Number of atoms to calculate = {}'.format(len(positions)))

    radial_distribution_values, distance_values = radial_distribution(distances, len(positions),
                                                                      minimum_distance=minimum_distance,
                                                                      maximum_distance=maximum_distance,
                                                                      verbose=verbose)
    pair_distribution_values = pair_distribution(distance_values, radial_distribution_values, mass_density)


    import matplotlib.pyplot as plt
    plt.xlabel(r'$r ({\rm \AA})$')
    plt.ylabel(r'$g(r) $')
    plt.xlim([minimum_distance, maximum_distance])
    text_x = minimum_distance + (maximum_distance - minimum_distance) * 0.05
    text_y = 0.9 * np.amax(pair_distribution_values)
    plt.text(text_x, text_y, label)

    if show_plot:
        plt.show()
    else:
        plt.savefig(filename)
        plt.close()
    return


def parse_command_line_arguments():
    """
    Take matdyn input file name from command line if provided
    :return: matdyn input filename as string
    """
    import argparse
    parser = argparse.ArgumentParser(description='Calculate PDF using pw.x output file positions...')
    parser.add_argument('-o', metavar="(pwscf.out)", type=str,
                        help='name of the pw.x output file', default='pwscf.out')
    parser.add_argument('--verbose', '-v', action=argparse.BooleanOptionalAction)
    parser.add_argument('--save_plot', '-s', action=argparse.BooleanOptionalAction)
    parser.add_argument('--maximum_distance', '-m', metavar="(maximum distance)", type=float, 
                        help='maximum distance [in angstroms] to plot g(r)', default=0.05)
    parser.add_argument('--density', '-d', metavar="(mass density)", type=float, 
                        help='reference mass density [in kg/m^3] to calculate g(r)', default=1000.)
    arguments = parser.parse_args()
    return arguments.o, arguments.verbose, arguments.save_plot, arguments.maximum_distance, arguments.density


def multislab_builder(positions,lx=3,ly=3):
    one_xarray = np.ones((len(positions),len(positions[0])))
    one_xarray[0:,1:]=0
    tmp_xarray = 1*positions
    empty_array = np.empty((0,3))
    empty_array=np.append(empty_array,tmp_xarray,axis=0)

    for i in np.arange(1,lx+1):
        tmp_xarray=tmp_xarray+one_xarray
        empty_array=np.append(empty_array,tmp_xarray,axis=0)

    one_yarray = np.ones((len(empty_array),len(empty_array[0])))
    one_yarray[0:,(0,-1)]=0
    tmp_yarray = 1*empty_array

    for i in np.arange(1,ly+1):
        tmp_yarray=tmp_yarray+one_yarray
        empty_array=np.append(empty_array,tmp_yarray,axis=0)
    return empty_array