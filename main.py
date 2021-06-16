from parameters import set_parameters
from data import get_data
from graph import graph_data


def main():
    is_slab, number_of_layers, exchange_correlation = set_parameters()
    data = get_data(number_of_layers, exchange_correlation)
    graph_data(data)

    return


if __name__ == '__main__':
    main()
