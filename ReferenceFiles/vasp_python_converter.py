def vasp_data_modifier(N_ML,phase,proj='n'):
    import numpy as np

    path = 'C:/Users/JOELA/PycharmProjects/silicon_nanomembrane_phase_transition/VestaFiles/'

    if phase == 'diamond':
        prefix = 'Si.Fd-3m_'
    if phase == 'betasn':
        prefix = 'Si.I4_1amd_'

    suffix = 'ML.POSCAR.VASP'


    with open(path+prefix+str(N_ML)+suffix) as f:
        lines = f.readlines()
    f.close()


    cella = float(lines[1])
    cellc = float(lines[4].split()[2])*cella
    nat = float(lines[6])


    position_lines = lines[8:int(8+nat)]
    positions = np.zeros((int(nat),3))

    i=0
    for x in position_lines:
        posq_string = x.split()
        posq = [float(q) for q in posq_string]
        pos = np.array([posq[0]*cella,posq[1]*cella,posq[2]*cellc])
        positions[i] = pos
        i += 1

    if proj == 'x':
        positions[:,1] = 0.0
        positions[:,2] = 0.0
    if proj == 'y':
        positions[:,0] = 0.0
        positions[:,2] = 0.0
    if proj == 'z':
        positions[:, 0] = 0.0
        positions[:, 1] = 0.0

    return positions, cella, cellc


def bulk_cell_vstacker(nml,phase,proj='n'):
    import numpy as np
    positions0 = vasp_data_modifier(0,phase,proj=proj)[0]
    cellc0 = vasp_data_modifier(0,phase,proj=proj)[2]

    positions= 1*positions0
    for n in np.arange(nml-1):
        column3_mod = positions0[:,2]+(n+1)*cellc0
        positionsn = 1*positions0
        positionsn[:,2] = column3_mod
        positions = np.concatenate((positions,positionsn))
    return positions


def ztest(n_atom,positions):
    import numpy as np

    zpos = positions[:, 2]

    zpos_sort_idxs = np.argsort(zpos)

    zpos_sort = zpos[zpos_sort_idxs]

    n_atom_zpos = zpos_sort[n_atom-1]

    return n_atom_zpos


def ind_atom_distances(n_atom,positions):
    import numpy as np

    check_for_array = (isinstance(n_atom, (list, tuple, np.ndarray)))
    if check_for_array == False:
        n_atom = np.array([n_atom])
    zpos = positions[:, int(2)]

    zpos_sort_idxs = np.argsort(zpos)

    positions_sort = positions[zpos_sort_idxs]

    distances = np.array([])
    for i in n_atom:
        atom_position = positions_sort[i-1]
        other_positions = np.delete(positions_sort,(i-1),axis=0)

        for j in other_positions:
            separation = atom_position - j
            distance = np.sqrt(np.sum(separation ** 2))
            distances = np.append(distances,distance)

    return np.sort(distances)

import numpy as np
