def add_displacements(coordinates, displacements):
    """
    Assume displacements have same units as coordinates
    Displacements and coordinates should be Numpy arrays of identical dimensions
    """
    import numpy as np

    if isinstance(displacements, np.ndarray) and isinstance(coordinates, np.ndarray):
        if displacements.shape == coordinates.shape:
            new_coordinates = displacements + coordinates
            return new_coordinates
        else:
            print('displacements and coordinates do not have same shape')
    else:
        print('displacements and coordinates are not both NumPy arrays')


if __name__ == '__main__':
    import numpy as np

    source_directory = '/Users/williamparker/Runs/QE/Si/diamond/pbe_uspp/phonons'
    source_scf_output = 'Si.Fd-3m.8-atom.scf.out'
    source_modes = 'matdyn.modes'

    from process_quantum_espresso_outputs import get_celldm, get_atom_types, get_positions, get_displacements

    source_atom_types = get_atom_types(source_directory + '/' + source_scf_output)
    source_coordinates = get_positions(source_directory + '/' + source_scf_output)
    lattice_parameter = get_celldm(source_directory + '/' + source_scf_output)
    phonon_displacements = get_displacements(source_directory + '/' + source_modes, q_point_index=2)

    ending_coordinates = add_displacements(source_coordinates, np.array(phonon_displacements))
    for atom_type, coordinate in zip(source_atom_types, ending_coordinates):
        print(' {} {:12.8f}{:12.8f}{:12.8f}'.format(atom_type, coordinate[0], coordinate[1], coordinate[2]))

