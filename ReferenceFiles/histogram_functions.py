def plot_histogram(layers,phase,minimum_distance=1.5,maximum_distance = 6,exchange_correlation = 'PBE',proj='n',individual_atom_plot=False,plot_n_atoms=1):
    import numpy as np
    from ReferenceFiles.vasp_python_converter import vasp_data_modifier, bulk_cell_vstacker
    from ReferenceFiles.pair_distribution import format_pdf,plot_pdf
    from ReferenceFiles.plot_formating import rgbcode_diamond,rgbcode_betasn,adjust_lightness,universal_linewidth,tick_fontsize,axis_fontsize
    import matplotlib.pyplot as plt

    if phase == 'diamond':
        prefix = 'Si.Fd-3m_'
        density = 2300
        rgbcode = rgbcode_diamond
    if phase == 'betasn':
        prefix = 'Si.I4_1amd_'
        density = 3370
        rgbcode = rgbcode_betasn


    n_entries = len(layers)
    k = 0


    for i in layers:   #np.array([0,1,2,h_n_ml]):
        output_file = 'DataFolder'+'/'+exchange_correlation+'/'+'output_files'+'/'+prefix+str(i)+'ML.scf.out'

        lattice_parameter = vasp_data_modifier(i,phase,proj=proj)[1]
        atomic_positions = vasp_data_modifier(i,phase,proj=proj)[0]

        if i == 0:
            atomic_positions = bulk_cell_vstacker(np.amax(layers),phase,proj=proj)


        if plot_n_atoms[0] == 'cut':
            nat_cut = plot_n_atoms[1]
            atoms = np.arange(len(atomic_positions))+1
            final_atoms = atoms[nat_cut:len(atoms)-nat_cut]
        elif plot_n_atoms[0] == 'rcut':
            nat_cut = plot_n_atoms[1]
            atoms = np.arange(len(atomic_positions))+1
            final_atoms = np.append(atoms[0:nat_cut],atoms[len(atoms)-nat_cut:len(atoms)])
        else:
            final_atoms = plot_n_atoms


        if individual_atom_plot == True:
            distance_valuesn, pair_distribution_valuesn = format_pdf(atomic_positions,density, minimum_distance=minimum_distance,maximum_distance=maximum_distance,individual_atom_histogram=True,n_atoms=final_atoms)
        else:
            distance_valuesn, pair_distribution_valuesn = format_pdf(atomic_positions,density, minimum_distance=minimum_distance,maximum_distance=maximum_distance)


        if maximum_distance <= 0.05:
            maximum_distance = np.ceil(lattice_parameter)


        if i == 0:
            rgbacode = adjust_lightness(rgbcode,0.5)
            plot_pdf(distance_valuesn, pair_distribution_valuesn,zorder=n_entries+1,color=rgbacode,linewidth=universal_linewidth,dashes=(3,3),label='Bulk')


        else:
            minalpha = 2.3
            topalpha = 0.7
            pdiv = (1/(n_entries))
            alpha = minalpha+k*pdiv*(topalpha-minalpha)
            rgbacode = adjust_lightness(rgbcode,alpha)
            plot_pdf(distance_valuesn, pair_distribution_valuesn,linestyle='solid',zorder=n_entries-k,color=rgbacode,linewidth=universal_linewidth,label=str(i)+'-ML')
        k+=1


    plt.xlabel(r'$r ({\rm \AA})$',fontsize=axis_fontsize)
    plt.ylabel(r'${g(r)}_{{N}_{ML}} $',fontsize=axis_fontsize)
    plt.xticks(fontsize=tick_fontsize)
    plt.yticks(fontsize=tick_fontsize)
    plt.xlim([minimum_distance, maximum_distance])
    return

def plot_histogram_differences(layers,phase,minimum_distance=1.5,maximum_distance = 6,exchange_correlation = 'PBE',proj='n',individual_atom_plot=False,plot_n_atoms=1):
    import numpy as np
    from ReferenceFiles.vasp_python_converter import vasp_data_modifier, bulk_cell_vstacker
    from ReferenceFiles.pair_distribution import format_pdf,plot_pdf
    from ReferenceFiles.plot_formating import rgbcode_diamond,rgbcode_betasn,adjust_lightness,universal_linewidth,tick_fontsize,axis_fontsize
    import matplotlib.pyplot as plt

    if phase == 'diamond':
        prefix = 'Si.Fd-3m_'
        density = 2300
        rgbcode = rgbcode_diamond
    if phase == 'betasn':
        prefix = 'Si.I4_1amd_'
        density = 3370
        rgbcode = rgbcode_betasn


    n_entries = len(layers)
    k = 0

    for i in layers:   #np.array([0,1,2,h_n_ml]):
        output_file = 'DataFolder'+'/'+exchange_correlation+'/'+'output_files'+'/'+prefix+str(i)+'ML.scf.out'
        lattice_parameter = vasp_data_modifier(i,phase,proj=proj)[1]
        atomic_positions = vasp_data_modifier(i,phase,proj=proj)[0]
        lattice_parameter0 = vasp_data_modifier(0,phase,proj=proj)[1]
        atomic_positions0 = bulk_cell_vstacker(i,phase,proj=proj)
        # atomic_positions0 = bulk_cell_vstacker(h_n_ml,phase,proj=proj)


        if i == 0:
            atomic_positions = bulk_cell_vstacker(1,phase,proj=proj)
            # atomic_positions = bulk_cell_vstacker(h_n_ml,phase,proj=proj)

        if plot_n_atoms[0] == 'cut':
            nat_cut = plot_n_atoms[1]
            atoms = np.arange(len(atomic_positions)) + 1
            final_atoms = atoms[nat_cut:len(atoms) - nat_cut]
        elif plot_n_atoms[0] == 'rcut':
            nat_cut = plot_n_atoms[1]
            atoms = np.arange(len(atomic_positions)) + 1
            final_atoms = np.append(atoms[0:nat_cut], atoms[len(atoms) - nat_cut:len(atoms)])
        else:
            final_atoms = plot_n_atoms


        if individual_atom_plot == True:
            distance_valuesn, pair_distribution_valuesn = format_pdf(atomic_positions,density, minimum_distance=minimum_distance,maximum_distance=maximum_distance,individual_atom_histogram=True,n_atoms=final_atoms)

            distance_values0, pair_distribution_values0 = format_pdf(atomic_positions0, density,minimum_distance=minimum_distance,maximum_distance=maximum_distance,individual_atom_histogram=True,n_atoms=final_atoms)

        else:
            distance_valuesn, pair_distribution_valuesn = format_pdf(atomic_positions,density, minimum_distance=minimum_distance,maximum_distance=maximum_distance)
            distance_values0, pair_distribution_values0 = format_pdf(atomic_positions0, density,minimum_distance=minimum_distance,maximum_distance=maximum_distance)


        pair_distribution_values = pair_distribution_valuesn - pair_distribution_values0


        if maximum_distance <= 0.05:
            maximum_distance = np.ceil(lattice_parameter)

        if i == 0:

            rgbacode = adjust_lightness(rgbcode,0.5)
            plot_pdf(distance_valuesn, pair_distribution_values0-pair_distribution_values0,zorder=n_entries+1,color=rgbacode,linewidth=universal_linewidth,dashes=(3,3),label='Bulk')

        else:
            minalpha = 2.3
            topalpha = 0.7
            pdiv = (1/(n_entries))
            alpha = minalpha+k*pdiv*(topalpha-minalpha)
            rgbacode = adjust_lightness(rgbcode,alpha)
            plot_pdf(distance_valuesn, pair_distribution_values,linestyle='solid',zorder=k+1,color=rgbacode,linewidth=universal_linewidth,label=str(i)+'-ML')
        k+=1


    plt.xlabel(r'$r ({\rm \AA})$',fontsize=axis_fontsize)
    plt.ylabel(r'${g(r)}_{{N}_{ML}}-{g(r)}_{Bulk} $',fontsize=axis_fontsize)
    plt.xticks(fontsize=tick_fontsize)
    plt.yticks(fontsize=tick_fontsize)
    plt.xlim([minimum_distance, maximum_distance])
    return

