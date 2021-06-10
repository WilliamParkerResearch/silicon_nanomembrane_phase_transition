
def parse_log_files():
    from DataFolder.PBE.relaxation.parse_log import parse_log_file
    from os import listdir
    from os.path import isfile, join

    # Location of log files
    directory = 'DataFolder/PBE/relaxation'
    directory_list = listdir(directory)

    # Get list of log files
    log_files = []
    for file in directory_list:
        if isfile(join(directory, file)) and file.split('.')[-1] == 'log':
            log_files.append(file)
    log_files.sort()

    total_energies, total_forces, phonon_frequencies = [], [], []
    for file in log_files:
        force, energy, frequencies = parse_log_file(directory + '/' + file)
        total_energies.append(energy)
        total_forces.append(force)
        phonon_frequencies.append(frequencies)

    return total_energies, total_forces, phonon_frequencies


def convert_to_electron_volt_angstrom_units(total_energies, total_forces, phonon_frequencies):
    import numpy as np

    from scipy.constants import value
    rydberg_to_electron_volt = value('Rydberg constant times hc in eV')
    bohr_to_angstrom = value('Bohr radius') * 10**10

    total_energies = np.array(total_energies) * rydberg_to_electron_volt
    total_forces = np.array(total_forces) * rydberg_to_electron_volt / bohr_to_angstrom
    phonon_frequencies = np.array(phonon_frequencies)

    return total_energies, total_forces, phonon_frequencies


def calculate_difference_from_converged_value(values):
    # values must by an NumPy ndarray
    import numpy as np

    if values.ndim == 2:
        # assuming nesting only one layer deep
        differences = []
        for internal_list in values:
            internal_differences = []
            for value in internal_list[:len(internal_list)-1]:
                internal_differences.append(value - internal_list[-1])
            differences.append(internal_differences)
    elif values.ndim == 1:
        differences = []
        for value in values[:len(values)-1]:
            differences.append(value - values[-1])
    else:
        print("Warning: Returning zeroes as differences")
        return np.zeros(values.shape)
    return np.array(differences)


def plot_energy_and_frequency_convergence_with_relaxation(total_energies, total_energy_differences, total_forces,
                                                          phonon_frequency_differences, figure_name=''):
    import matplotlib.pyplot as plt
    figure, axes = plt.subplots(nrows=2, ncols=1, sharex=True, gridspec_kw={'height_ratios': [1, 2]})
    # Delete space between plots
    figure.subplots_adjust(hspace=0)

    # Energy convergence plot
    axes[0].scatter(total_forces[:len(total_energies)-1], total_energy_differences, marker='o')
    axes[0].invert_xaxis()
    #axes[0].set_xlabel(r'$F_{\rm total}$ (eV/Å)')
    axes[0].set_ylabel(r'$\Delta E_{\rm total}$ (eV)')
    axes[0].set_yscale('log')

    # Frequency convergence plot
    for mode_differences in phonon_frequency_differences:
        axes[1].scatter(total_forces[:len(total_energies)-1], mode_differences, marker='o')

    axes[1].set_ylabel(r'$\Delta \omega$ (THz)')
    axes[1].set_ylim([-1.e0, 1.e0])
    axes[1].set_yscale('symlog', linthresh=.01)
    axes[1].set_yticks([-0.1, -0.01, 0, 0.01, 0.1])
    axes[1].axhline(color='black')

    axes[1].invert_xaxis()
    axes[1].set_xlabel(r'$F_{\rm total}$ (eV/Å)')

    for axis in axes:
        axis.set_xscale('log')
        axis.set_xlim([1.e0, 1.e-3])

    if len(figure_name) > 0:
        plt.savefig(figure_name)

    plt.show()

    return


if __name__ == '__main__':
    energies, forces, frequencies = parse_log_files()
    energies, forces, frequencies = convert_to_electron_volt_angstrom_units(energies, forces, frequencies)

    energy_differences = calculate_difference_from_converged_value(energies)
    frequency_differences = calculate_difference_from_converged_value(frequencies.T)

    figure_name = 'Si.PBE.Fd-3m_1ML.energy_frequency_relaxation.png'
    plot_energy_and_frequency_convergence_with_relaxation(energies, energy_differences, forces, frequency_differences,
                                                          figure_name=figure_name)
