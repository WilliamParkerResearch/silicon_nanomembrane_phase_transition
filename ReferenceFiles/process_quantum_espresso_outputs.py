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


if __name__ == '__main__':
    filename = '/Users/williamparker/Runs/QE/Si/1_ML/pbe_uspp/edos/Si.1_ML.scf.out'
    lattice_parameter = get_celldm(filename)
    print('Found celldm(1) = {}'.format(lattice_parameter))
    atomic_positions = lattice_parameter * get_positions(filename)
    print('Found taus:')
    for position in atomic_positions:
        print(position)
