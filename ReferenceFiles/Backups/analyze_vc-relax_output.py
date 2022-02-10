import numpy as np


def get_values(output_file, search_key, key_placement_index, value_placement_index, number_of_values=5):
    """
        Parameters
        ----------
           output_file : str
                SCF vc-relax output file name and path
           search_key : str
                name to search for in output file
           key_placement_index : int
                index in line split to look for key
           value_placement_index : int
                index in line split to get value from
           number_of_values : int, optional
                number of values from output file to report on
        Returns
        -------
           values : ndarray(number_of_values)
                values extracted from output file
    """
    values = []
    with open(output_file) as out_file:
        for line in out_file.readlines():
            if len(line.split()) >= value_placement_index:
                if line.split()[key_placement_index] == search_key:
                    values.append(float(line.split()[value_placement_index]))
    if number_of_values > 0:
        report_index = len(values) - number_of_values
    elif number_of_values == 0:
        report_index = 0
    return np.array(values[report_index:])


def print_table_of_values(output_file, number_of_values_to_print=4):
    total_energies = get_values(output_file, '!', 0, 4, number_of_values=number_of_values_to_print+1)
    total_forces = get_values(output_file, 'Total', 0, 3, number_of_values=number_of_values_to_print+1)
    cell_pressures = get_values(output_file, 'P=', 4, 5, number_of_values=number_of_values_to_print+1)

    change_in_energy = np.diff(total_energies)
    change_in_total_force = np.diff(total_forces)
    change_in_pressure = np.diff(cell_pressures)

    dashed_line = 120 * '-'
    print('{}'.format(dashed_line))

    energy_string = 'E_total (Ry)'
    delta_energy_string = 'ΔE_total (Ry)'
    force_string = 'F_total (Ry/bohr)'
    delta_force_string = 'ΔF_total (Ry/bohr)'
    pressure_string = 'P (kbar)'
    delta_pressure_string = 'ΔP (kbar)'
    print('{:20s}{:20s}{:20s}{:20s}{:20s}{:20s}'.format(energy_string, delta_energy_string, force_string,
                                                        delta_force_string, pressure_string, delta_pressure_string))
    print('{}'.format(dashed_line))

    for energy, delta_e, force, delta_f, pressure, delta_p in \
            zip(total_energies[1:], change_in_energy, total_forces[1:], change_in_total_force,
                cell_pressures[1:], change_in_pressure):
        print('{:12.8f}  {:18.8f}  {:18.6f}  {:18.6f}  {:12.2f}  {:18.2f}\n'.format(
            energy, delta_e, force, delta_f, pressure, delta_p))

    print('{}'.format(dashed_line))
    return


def graph_values(output_file, number_of_values_to_graph=-1):
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(3)

    total_energies = get_values(output_file, '!', 0, 4, number_of_values=number_of_values_to_graph+1)
    indices = np.arange(1, len(total_energies)+1)
    axes[0].plot(indices, total_energies)
    axes[0].set(ylabel=r'$E_{total}$ (Ry)')

    total_forces = get_values(output_file, 'Total', 0, 3, number_of_values=number_of_values_to_graph+1)
    axes[1].plot(indices, total_forces, color='orange')
    axes[1].set(ylabel=r'$F_{total}$ (Ry/bohr)')

    cell_pressures = get_values(output_file, 'P=', 4, 5, number_of_values=number_of_values_to_graph+1)
    axes[2].plot(indices, cell_pressures, color='green')
    axes[2].set(xlabel='Relaxation step', ylabel=r'$P$ (kbar)')
    axes[2].axhline(color='black', alpha=0.5, linestyle='dashed')

    from matplotlib.ticker import FormatStrFormatter
    for axis in axes.flat:
        axis.label_outer()
        axis.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))

    figure.subplots_adjust(left=0.2, hspace=0.)
    plt.show()
    return


if __name__ == '__main__':
    filename = '/Users/williamparker/Documents/Research/Silicon Nanomembrane Phase Transition/output_files/Si.1_ML.Fd-3m.vc-relax.out.8'
    print_table_of_values(filename)
    graph_values(filename)

