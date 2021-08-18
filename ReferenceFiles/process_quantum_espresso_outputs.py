import numpy as np


def get_celldm(output_file, output_type='scf', verbose=False):
    if verbose:
        print('Getting celldm from {} ...'.format(output_file))
        print('\t Output file {} is type {}'.format(output_file, output_type))
    if output_type == 'scf':
        with open(output_file) as out_file:
            for line in out_file.readlines():
                if len(line.split()) > 1:
                    if line.split()[0] == 'celldm(1)=':
                        celldm = float(line.split()[1])
    try:
        if verbose:
            print('\t Celldm {} bohr found'.format(celldm))
            print()
        return celldm
    except NameError:
        print('celldm(1) not found in {}'.format(output_type))


def get_positions(output_file, output_type='scf', verbose=False):
    """Returns atomic positions in alat units"""
    if verbose:
        print('Getting atomic positions from {} ...'.format(output_file))
        print('\t Output file {} is type {}'.format(output_file, output_type))
    positions = []
    if output_type == 'scf':
        with open(output_file) as out_file:
            for line in out_file.readlines():
                if len(line.split()) > 2:
                    if line.split()[2] == 'tau(':
                        position = [float(line.split()[6]), float(line.split()[7]), float(line.split()[8])]
                        positions.append(position)
    try:
        if verbose:
            print('\t Positions for {} atoms found'.format(len(positions)))
            print('\t Position 1 [alat] = {} '.format(positions[0]))
            print()
        return np.array(positions)
    except NameError:
        print(' Positions lines (tau) not found in {}'.format(output_type))


def get_atom_types(output_file, output_type='scf', verbose=False):
    """Returns array of atom types"""
    if verbose:
        print('Getting atom types from {} ...'.format(output_file))
        print('\t Output file {} is type {}'.format(output_file, output_type))
    atom_types = []
    if output_type == 'scf':
        with open(output_file) as out_file:
            for line in out_file.readlines():
                if len(line.split()) > 2:
                    if line.split()[2] == 'tau(':
                        atom_type = line.split()[1]
                        atom_types.append(atom_type)
    try:
        if verbose:
            print('\t Atom types for {} atoms found'.format(len(atom_types)))
            print('\t Atom type 1 = {} '.format(atom_types[0]))
            print()
        return atom_types
    except NameError:
        print(' Positions lines (tau) not found in {}'.format(output_type))


def count_qpoints(modes_file):
    # Count the number of q-points in a matdyn-generated modes file
    number_of_q_points = 0
    with open(modes_file) as file:
        for line in file.readlines():
            if len(line.split()) > 0 and line.split()[0] == 'q':
                number_of_q_points += 1
    return number_of_q_points


def count_modes(modes_file):
    number_of_modes = 0
    frequency_index = 1
    with open(modes_file) as file:
        for line in file.readlines():
            if len(line.split()) > 0 and line.split()[0] == 'freq' and line.split()[2] == str(frequency_index) + ')':
                frequency_index += 1
                number_of_modes += 1
    return number_of_modes


def get_displacements(modes_file, q_point_index=1, frequency_index=1):
    # q-points and frequencies indexed starting at 1, but
    # q-point and frequency line number indices start at 0
    number_of_modes = count_modes(modes_file)
    number_of_atoms = int(number_of_modes/3)   # each atom has three spatial degrees of freedom
    print('Expecting {} modes for {} atoms'.format(number_of_modes, number_of_atoms))
    q_point_line_number = 2 + (q_point_index - 1) * ((number_of_atoms + 1) * number_of_modes + 5)
    frequency_line_number = q_point_line_number + (frequency_index-1) * (number_of_atoms + 1) + 2
    print('Looking for frequency number {} in q-point number {}'.format(frequency_index, q_point_index))
    print(' at line {}, located after the q-point in line {}'.format(frequency_line_number, q_point_line_number))

    with open(modes_file) as file:
        displacements = []
        q_point_found = False
        frequency_found = False
        for line_index, line in enumerate(file.readlines()):
            if len(line.split()) > 0:
                if line_index == q_point_line_number and line.split()[0] == 'q':
                    q_point = [ float(line.split()[2]), float(line.split()[3]), float(line.split()[4]) ]
                    print('Q-point # {} = {}'.format(q_point_index, q_point))
                    q_point_found = True
                    continue
                if q_point_found:
                    if line_index == frequency_line_number and line.split()[0] == 'freq':
                        print('\t Frequency {} found in line {}'.format(frequency_index, line_index))
                        print('\t Adding displacements for ω({}) = {} THz'.format(frequency_index, line.split()[4]))
                        print()
                        frequency_found = True
                        continue
                    if frequency_found:
                        if frequency_line_number < line_index <= (frequency_line_number + number_of_atoms)\
                                and line.split()[0] == '(':
                            displacement = [float(line.split()[1]), float(line.split()[3]), float(line.split()[5])]
                            displacements.append(displacement)
    return np.array(displacements)


if __name__ == '__main__':
    filename = '/Users/williamparker/PycharmProjects/silicon_nanomembrane_phase_transition/DataFolder/PBE/modes/matdyn.modes'
    atomic_displacements = get_displacements(filename, q_point_index=301, frequency_index=24)
    print(atomic_displacements)

