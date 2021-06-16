def set_parameters():
    # Parameters
    is_slab = False
    if is_slab:
        number_of_layers = 1
    else:
        number_of_layers = 0
    exchange_correlation = 'PBE'

    # Print parameters
    print('Calculation parameters:')
    print('\t Exchange-correlation = {}'.format(exchange_correlation))
    if is_slab:
        print('\t {}-layer slabs')
    else:
        print('\t Bulk structures')

    return is_slab, number_of_layers, exchange_correlation
