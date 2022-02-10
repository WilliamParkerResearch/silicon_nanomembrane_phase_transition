def pdos(fname,type='r'):    # use r for individual atoms and p for tot
    import numpy as np

    data = np.loadtxt(fname)

    k = np.unique(data[:,0])
    e = np.unique(data[:,1])
    dos = np.zeros([len(k), len(e)])

    for i in range(len(data)):
        e_index = int(i % len(e))
        k_index = int(data[i][0] - 1)
        if type == 'r':
            dos[k_index, e_index] = data[i][2]
        elif type == 'p':
            dos[k_index, e_index] = data[i][3]

    pdos = np.zeros([1,len(dos[0])])[0]

    for i in range(len(dos)):
        pdos= pdos+dos[i]
    return e,pdos