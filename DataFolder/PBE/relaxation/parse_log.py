import numpy as np


def parse_log_file(filename):
    with open(filename) as log_file:
        first_line = log_file.readline()
        total_force = float(first_line.split()[4])
        second_line = log_file.readline()
        total_energy = float(second_line.split()[4])
        frequencies = []
        for line in log_file.readlines():
            frequencies.append(float(line.split()[4]))
    return total_force, total_energy, frequencies


if __name__ == '__main__':
    force, energy, phonon_frequencies = parse_log_file('max1.log')
