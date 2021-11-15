def eos_properties(N_ML,number_of_volume_points=10000, n_pressure_values=1000000, vol_ext=0.5, exchange_correlation='PBE', OF=False, etype='ce', n_etype=0, eos_type = 'murnaghan'):
    import matplotlib.pyplot as plt
    import scipy.optimize as sp
    from scipy import interpolate
    import numpy as np
    from importlib import import_module
    from ReferenceFiles.enthalpy import enthalpy_from_volume
    from ReferenceFiles.FunctionDefinitions import mid,square_differences
    from ReferenceFiles.equations_of_state import birch_murnaghan,murnaghan,vinet,pressure_from_energy_equation_of_state

    if OF == False:
        directoryofdata = 'DataFolder' + '.' + exchange_correlation + '.' + 'eostm' + '.' + 'Data_' + str(N_ML) + 'L'
    else:
        directoryofdata = 'DataFolder' + '.' + exchange_correlation + '.' + 'eost' + '.' + 'Data_' + str(N_ML) + 'L' + '_' + etype + str(n_etype)
    data = import_module(directoryofdata)

    if eos_type == 'birch-murnaghan':
        def eos(p, v):
            return birch_murnaghan(p, v)
    elif eos_type == 'murnaghan':
        def eos(p, v):
            return murnaghan(p, v)
    elif eos_type == 'vinet':
        def eos(p, v):
            return vinet(p, v)


    def intersection_idx(equation):  # prints lowerbound index in which the value is found within it and the index above it
        i = np.argwhere(np.diff(np.sign(equation))).flatten()
        return i


    def extend_array_extremas(array, percentage=0):
        min_value = np.amin(array)
        max_value = np.amax(array)
        value_diff = (percentage) * np.absolute(max_value - min_value)
        return value_diff


    #   Diamond structure calculations
    volumes_sim_diamond = data.volumes_sim_diamond
    total_energies_strain_diamond = data.total_energies_strain_diamond
    initial_parameters_diamond = (total_energies_strain_diamond[mid(total_energies_strain_diamond)], 1e11, 3.5,
                                  volumes_sim_diamond[mid(volumes_sim_diamond)])
    fit_parameters_diamond = sp.fmin(square_differences, initial_parameters_diamond,
                                     args=(volumes_sim_diamond, total_energies_strain_diamond, eos), maxiter=100000)
    volume_diamond_extension = extend_array_extremas(volumes_sim_diamond, vol_ext)
    volumes_diamond = np.linspace(np.amin(volumes_sim_diamond) - volume_diamond_extension,
                                  np.amax(volumes_sim_diamond) + volume_diamond_extension, num=number_of_volume_points)
    pressures_diamond = -pressure_from_energy_equation_of_state(fit_parameters_diamond,volumes_diamond,eos=eos_type)
    enthalpy_diamond = enthalpy_from_volume(fit_parameters_diamond, volumes_diamond, eos=eos_type)
    fit_total_energies_strain_diamond = eos(fit_parameters_diamond, volumes_diamond)
    #     Beta-Sn structure calculations
    volumes_sim_betasn = data.volumes_sim_betasn
    total_energies_strain_betasn = data.total_energies_strain_betasn
    initial_parameters_shape_betasn = (total_energies_strain_betasn[mid(total_energies_strain_betasn)], 1e11, 3.5,
                                       volumes_sim_betasn[mid(volumes_sim_betasn)])
    fit_parameters_betasn = sp.fmin(square_differences, initial_parameters_shape_betasn,
                                    args=(volumes_sim_betasn, total_energies_strain_betasn, eos), maxiter=100000)
    volume_betasn_extension = extend_array_extremas(volumes_sim_betasn, vol_ext)
    volumes_betasn = np.linspace(np.amin(volumes_sim_betasn) - volume_betasn_extension,
                                 np.amax(volumes_sim_betasn) + volume_betasn_extension, num=number_of_volume_points)
    pressures_betasn = -pressure_from_energy_equation_of_state(fit_parameters_betasn,volumes_betasn,eos=eos_type)
    enthalpy_betasn = enthalpy_from_volume(fit_parameters_betasn, volumes_betasn, eos=eos_type)
    fit_total_energies_strain_betasn = eos(fit_parameters_betasn, volumes_betasn)
    #           parameters = (E0, K0, K0', V0)

    idx_d = np.argsort(pressures_diamond)
    pressures_diamond = pressures_diamond[idx_d]
    enthalpy_diamond = enthalpy_diamond[idx_d]

    idx_b = np.argsort(pressures_betasn)
    pressures_betasn = pressures_betasn[idx_b]
    enthalpy_betasn = enthalpy_betasn[idx_b]


    pressure_min = np.amax(np.array([pressures_diamond[0],pressures_betasn[0]]))
    pressure_max = np.amin(np.array([pressures_diamond[-1],pressures_betasn[-1]]))
    pressures = np.linspace(pressure_min,pressure_max,n_pressure_values)

    enthalpy_diamond_func = interpolate.interp1d(pressures_diamond,enthalpy_diamond)
    enthalpy_diamond_fit = enthalpy_diamond_func(pressures)
    enthalpy_betasn_func = interpolate.interp1d(pressures_betasn,enthalpy_betasn)
    enthalpy_betasn_fit = enthalpy_betasn_func(pressures)

    enthalpy_diff = enthalpy_diamond_fit-enthalpy_betasn_fit

    p_idx = intersection_idx(enthalpy_diff)
    t_pressure = pressures[int(p_idx)]

    t_vol_d_idx = intersection_idx(pressures_diamond-t_pressure)
    t_volume_diamond = float(volumes_diamond[t_vol_d_idx])

    t_vol_b_idx = intersection_idx(pressures_betasn-t_pressure)
    t_volume_betasn = float(volumes_betasn[t_vol_b_idx])
    return t_pressure,t_volume_diamond,t_volume_betasn,fit_parameters_diamond,fit_parameters_betasn


