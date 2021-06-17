import matplotlib.pyplot as plt
# import numpy as np
from ReferenceFiles.unit_conversions import cubic_meters_per_cubic_angstrom, eV_per_joule


def graph_data(data):

    fig, ax = plt.subplots(ncols=3)

    # E vs V
    for structure in data.keys():
        volumes = data[structure]['volumes'] / cubic_meters_per_cubic_angstrom
        energies = data[structure]['energies'] * eV_per_joule
        ax[0].plot(volumes,
                 energies,
                 marker='o', linestyle='')

    # Fit E(V)
    ax[0].set_xlabel(r'$V$ (Å$^3$/atom)')
    ax[0].set_ylabel(r'$E$ (eV/atom)')

    # P(V)
    ax[1].set_xlabel(r'$V$ (Å$^3$/atom)')
    ax[1].set_ylabel(r'$P$ (GPa)')

    # H(V)
    ax[2].set_xlabel(r'$V$ (Å$^3$/atom)')
    ax[2].set_ylabel(r'$H$ (eV/atom)')

    plt.tight_layout()
    plt.savefig('test.png')

    return