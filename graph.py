import matplotlib.pyplot as plt
# import numpy as np


def graph_data(data):

    for structure in data.keys():
        plt.plot(data[structure]['volumes'],
                 data[structure]['energies'],
                 marker='o')

    plt.savefig('test.png')

    return