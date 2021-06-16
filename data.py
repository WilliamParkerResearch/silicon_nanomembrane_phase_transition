import importlib


def get_data(number_of_layers, exchange_correlation):
    # Get data for given parameters
    data_directory = 'DataFolder.' + str(exchange_correlation) + '.eos.'
    data_file = 'Data_' + str(number_of_layers) + 'L'
    data_module = importlib.import_module(data_directory + data_file)
    diamond_energies = data_module.total_energies_strain_diamond
    diamond_volumes = data_module.volumes_sim_diamond
    betasn_energies = data_module.total_energies_strain_betasn
    betasn_volumes = data_module.volumes_sim_betasn

    # Structure data for output
    data = {
        'diamond': {'energies': [], 'volumes': []},
        'betasn': {'energies': [], 'volumes': []}
    }
    data['diamond']['energies'] = diamond_energies
    data['diamond']['volumes'] = diamond_volumes
    data['betasn']['energies'] = betasn_energies
    data['betasn']['volumes'] = betasn_volumes

    return data


