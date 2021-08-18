def get_positions(filename):
    import re
    elements = []
    coordinates = []
    relax_flags = []
    with open(filename) as file:
        for line in file.readlines():
            line_length = len(line.split())
            if line_length == 2:
                if line.split()[0] == 'ATOMIC_POSITIONS':
                    coordinate_system = line.split()[1]
            if line_length == 4 or line_length == 7:
                match = re.search('[A-Z]', line.split()[0])
                if match:
                    elements.append(line.split()[0])
                    coordinate = [float(line.split()[1]), float(line.split()[2]), float(line.split()[3])]
                    coordinates.append(coordinate)
                    if line_length == 7:
                        relax_flag_set = [int(line.split()[4]), int(line.split()[5]), int(line.split()[6])]
                        relax_flags.append(relax_flag_set)

    return elements, coordinates, relax_flags, coordinate_system


def get_lattice_parameter(filename):
    with open(filename) as file:
        for line in file.readlines():
            if line.split()[0] == 'celldm(1)':
                lattice_parameter = line.split()[2]

    return lattice_parameter


if __name__ == '__main__':
    import numpy as np
    input_file_directory = '/Users/williamparker/Runs/QE/Si/1_ML/pbe_uspp/edos'
    input_file_name = 'Si.1_ML.scf.in'
    input_file = input_file_directory + '/' + input_file_name
    alat = get_lattice_parameter(input_file)
    print(alat)
    input_elements, input_coordinates, coordinate_relax_flags, input_coordinate_system = get_positions(input_file)
    print(input_elements)
    print(input_coordinates)
    print(input_coordinate_system)
