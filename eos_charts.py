from matplotlib.ticker import StrMethodFormatter
import numpy as np
# To Do
# Determine why transition pressures are wrong

# Parameters
system = 'slab'
# Import correct file
if system == 'slab':
    from eos_take2 import *
elif system == 'bulk':
    import matplotlib.pyplot as plt
    from eos_information import *
# exchange_correlation = 'PBE'          ###CURRENTLY NEED TO CHANGE IN eos_take2.py or eos_information.py
write_png_file = False                 # Write PNG file if True, open plot window if False
tangent_line_shift = -7.0e-2            # Shift tangent line by this factor of the vertical axis range
text_vertical_shift_factor = 5.0e-2     # Shift texts by this factor of the vertical axis range
text_horizontal_shift_factor = 5.0e-3   # Shift texts by this factor of the horizontal axis range
#   Color choices
diamond_color = 'blue'
beta_tin_color = 'red'
tangent_line_color = 'green'
#   Formatting
vertical_axis_tick_format = '{x:.1f}'
string_format = '{0:3.1f}'
#   Label-specific shift factors
if system == 'slab':
    if exchange_correlation == 'PZ':
        diamond_equilibrium_volume_horizontal_shift_factor = 10
        diamond_bulk_modulus_horizontal_shift_factor = 4
        beta_tin_bulk_modulus_horizontal_shift_factor = 0
        transition_pressure_vertical_shift_factor = 2
        transition_pressure_horizontal_shift_factor = 5
    elif exchange_correlation == 'PBE':
        diamond_equilibrium_volume_horizontal_shift_factor = 1
        diamond_bulk_modulus_horizontal_shift_factor = 4
        beta_tin_bulk_modulus_horizontal_shift_factor = 0
        transition_pressure_vertical_shift_factor = 2
        transition_pressure_horizontal_shift_factor = 15
    elif exchange_correlation == 'SCAN':
        diamond_equilibrium_volume_horizontal_shift_factor = 1
        diamond_bulk_modulus_horizontal_shift_factor = 6
        beta_tin_bulk_modulus_horizontal_shift_factor = 0
        transition_pressure_vertical_shift_factor = 2
        transition_pressure_horizontal_shift_factor = 10
elif system == 'bulk':
    if exchange_correlation == 'PZ':
        diamond_equilibrium_volume_horizontal_shift_factor = 10
        diamond_bulk_modulus_horizontal_shift_factor = 14
        beta_tin_bulk_modulus_horizontal_shift_factor = 22
        transition_pressure_vertical_shift_factor = 2
        transition_pressure_horizontal_shift_factor = 10
    elif exchange_correlation == 'PBE':
        diamond_equilibrium_volume_horizontal_shift_factor = 1
        diamond_bulk_modulus_horizontal_shift_factor = 14
        beta_tin_bulk_modulus_horizontal_shift_factor = 22
        transition_pressure_vertical_shift_factor = 1
        transition_pressure_horizontal_shift_factor = 5
    elif exchange_correlation == 'SCAN':
        diamond_equilibrium_volume_horizontal_shift_factor = 10
        diamond_bulk_modulus_horizontal_shift_factor = 18
        beta_tin_bulk_modulus_horizontal_shift_factor = 22
        transition_pressure_vertical_shift_factor = 2
        transition_pressure_horizontal_shift_factor = 5

#   Volume range
total_volumes=np.append(volumes_sim_BetaSn,volumes_sim_diamond)
minvol=int(np.amin(total_volumes)/cubic_meters_per_cubic_angstrom)
maxvol=int(np.amax(total_volumes)/cubic_meters_per_cubic_angstrom)
voldisp=(maxvol-minvol)*0.05
volume_range = (minvol-voldisp, maxvol+voldisp)
delta_volume = 5

# Use LaTeX
mpl.rcParams['text.usetex'] = True

# Font settings for presentation slide graphics
mpl.rcParams['font.serif'] = "Times"
mpl.rcParams['font.family'] = "Georgia"
mpl.rcParams['font.size'] = 18

# Set up plot
fig, ax = plt.subplots()

# Text strings to overlay on plot
equilibrium_volume_string = r'$V_0$ = '
volume_unit_string = r' \r{A}$^3$'
bulk_modulus_string = r'$K_0$ = '
transition_pressure_string = r'$P_T$ = '
pressure_unit_string = ' GPa'

# Unit conversions
#    Volumes
volumes_diamond = volumes_diamond/cubic_meters_per_cubic_angstrom
volumes_sim_diamond = volumes_sim_diamond/cubic_meters_per_cubic_angstrom
volumes_BetaSn = volumes_BetaSn/cubic_meters_per_cubic_angstrom
volumes_sim_BetaSn = volumes_sim_BetaSn/cubic_meters_per_cubic_angstrom
#    Energies
energy_conversion_factor = 6.242e18
fit_total_energies_strain_diamond = fit_total_energies_strain_diamond * energy_conversion_factor
total_energies_strain_diamond = total_energies_strain_diamond * energy_conversion_factor

if system == 'slab':
    fit_total_energies_strain_shape_BetaSn = fit_total_energies_strain_BetaSn * energy_conversion_factor
    total_energies_strain_shape_BetaSn = total_energies_strain_BetaSn * energy_conversion_factor
elif system == 'bulk':
    fit_total_energies_strain_shape_BetaSn = fit_total_energies_strain_shape_BetaSn * energy_conversion_factor
    total_energies_strain_shape_BetaSn = total_energies_strain_shape_BetaSn * energy_conversion_factor

phase_energy_difference = np.abs(np.amin(total_energies_strain_diamond) - np.amin(total_energies_strain_shape_BetaSn))

# Plots
#   Diamond
ax.plot(volumes_diamond, fit_total_energies_strain_diamond, color=diamond_color)
ax.scatter(volumes_sim_diamond, total_energies_strain_diamond, color=diamond_color)
#   Beta-tin
ax.plot(volumes_BetaSn, fit_total_energies_strain_shape_BetaSn, color=beta_tin_color)
ax.scatter(volumes_sim_BetaSn, total_energies_strain_shape_BetaSn, color=beta_tin_color, marker='s')
#   Tangent line
tangent_line_values = (murnaghan(fit_parameters_shape_BetaSn, tvol_beta)-(volume-tvol_beta)
                       * tpressure) * energy_conversion_factor
volume = volume/cubic_meters_per_cubic_angstrom
# diamond_minimum_difference_index = np.argmin(np.abs(tangent_line_values-fit_total_energies_strain_diamond))
# beta_tin_minimum_difference_index = np.argmin(np.abs(tangent_line_values-fit_total_energies_strain_shape_BetaSn))
# if beta_tin_minimum_difference_index < diamond_minimum_difference_index:
#     ax.plot(volume[:diamond_minimum_difference_index],
#             tangent_line_values[:diamond_minimum_difference_index] +
#             tangent_line_shift,
#             color=tangent_line_color)
# else:
#     ax.plot(volume[diamond_minimum_difference_index:beta_tin_minimum_difference_index],
#             tangent_line_values[diamond_minimum_difference_index:beta_tin_minimum_difference_index] +
#             tangent_line_shift,
#             color=tangent_line_color)
ax.plot(volume,
        tangent_line_values +
        tangent_line_shift,
        color=tangent_line_color)

# Axes properties
#   Horizontal
ax.set_xlim(left=volume_range[0], right=volume_range[1])
ax.set_xticks(np.arange(int(volume_range[0]), int(volume_range[1]), delta_volume))
ax.set_xlabel(r'Volume (\r{A}$^3$/atom)')
#   Vertical
ax.set_ylim(bottom=np.amin(total_energies_strain_diamond) - 1.5*phase_energy_difference,
            top=np.amax(total_energies_strain_shape_BetaSn) + 0.5*phase_energy_difference)
ax.set_yscale(r'linear')
ax.set_ylabel(r'Total Energy (eV/atom)')
ax.yaxis.set_major_formatter(StrMethodFormatter(vertical_axis_tick_format))

# Equilibrium volume line segments
diamond_equilibrium_volume = fit_parameters_diamond[3] / cubic_meters_per_cubic_angstrom
diamond_equilibrium_energy = fit_parameters_diamond[0] * energy_conversion_factor
beta_tin_equilibrium_volume = fit_parameters_shape_BetaSn[3] / cubic_meters_per_cubic_angstrom
beta_tin_equilibrium_energy = fit_parameters_shape_BetaSn[0] * energy_conversion_factor
plot_energy_minimum, plot_energy_maximum = ax.get_ylim()
plot_energy_range = plot_energy_maximum - plot_energy_minimum
diamond_line_segment_fraction = (diamond_equilibrium_energy-plot_energy_minimum)/plot_energy_range
beta_tin_line_segment_fraction = (beta_tin_equilibrium_energy-plot_energy_minimum)/plot_energy_range
ax.axvline(diamond_equilibrium_volume, ymax=diamond_line_segment_fraction, color='black', linestyle='--')
ax.axvline(beta_tin_equilibrium_volume, ymax=diamond_line_segment_fraction, color='black', linestyle='--')
# ax.axvline(beta_tin_equilibrium_volume, ymax=beta_tin_line_segment_fraction, color='black', linestyle='--')

# Text annotations
#    Locations
plot_volume_minimum, plot_volume_maximum = ax.get_xlim()
plot_volume_range = plot_volume_maximum - plot_volume_minimum

text_vertical_shift = text_vertical_shift_factor * plot_energy_range
text_horizontal_shift = text_horizontal_shift_factor * plot_volume_range

diamond_equilibrium_volume_coordinates = (diamond_equilibrium_volume +
                                          text_horizontal_shift * diamond_equilibrium_volume_horizontal_shift_factor,
                                          plot_energy_minimum + text_vertical_shift)
#       Vertically aligned with average first and last total energy
diamond_bulk_modulus_coordinates = (diamond_equilibrium_volume +
                                    text_horizontal_shift*diamond_bulk_modulus_horizontal_shift_factor,
                                    0.5*(total_energies_strain_diamond[0] +
                                         total_energies_strain_diamond[-1]))
#       Aligned with largest total energy
# diamond_bulk_modulus_coordinates = (diamond_equilibrium_volume,
#                                     np.amax(total_energies_strain_diamond))
beta_tin_equilibrium_volume_coordinates = (beta_tin_equilibrium_volume + text_horizontal_shift,
                                          plot_energy_minimum + text_vertical_shift)
#        Vertically aligned with average first and last total energy
beta_tin_bulk_modulus_coordinates = (beta_tin_equilibrium_volume +
                                     text_horizontal_shift*beta_tin_bulk_modulus_horizontal_shift_factor,
                                     0.5*(total_energies_strain_shape_BetaSn[0] +
                                          total_energies_strain_shape_BetaSn[-1]))
#        Aligned with largest total energy
# beta_tin_bulk_modulus_coordinates = (beta_tin_equilibrium_volume,
#                                      np.amax(total_energies_strain_shape_BetaSn))

if t_pres > 0:
    transition_pressure_coordinates = (plot_volume_minimum +
                                       transition_pressure_horizontal_shift_factor*text_horizontal_shift,
                                       diamond_equilibrium_energy -
                                       transition_pressure_vertical_shift_factor*text_vertical_shift)
    ax.text(transition_pressure_coordinates[0], diamond_equilibrium_volume_coordinates[1] - 0.25 * text_vertical_shift,
             r'$V_0$ = ', color='black')
else:
    transition_pressure_coordinates = (beta_tin_equilibrium_volume +
                                       transition_pressure_horizontal_shift_factor*text_horizontal_shift,
                                       total_energies_strain_diamond[-1] -
                                       transition_pressure_vertical_shift_factor*text_vertical_shift)
    ax.text(volumes_sim_diamond[0], diamond_equilibrium_volume_coordinates[1] - 0.25*text_vertical_shift,
             r'$V_0$ = ', color='black')


#   Strings
diamond_equilibrium_volume_string = string_format.format(diamond_equilibrium_volume) + \
                                    volume_unit_string
diamond_bulk_modulus_string = bulk_modulus_string + string_format.format(k_0_diamond) + \
                              pressure_unit_string
beta_tin_equilibrium_volume_string = string_format.format(beta_tin_equilibrium_volume) + \
                                    volume_unit_string
beta_tin_bulk_modulus_string = bulk_modulus_string + string_format.format(k_0_betasn) + \
                              pressure_unit_string
transition_pressure_string = transition_pressure_string + string_format.format(t_pres) + pressure_unit_string


#    Placements
ax.text(diamond_equilibrium_volume_coordinates[0], diamond_equilibrium_volume_coordinates[1],
        diamond_equilibrium_volume_string, color='black')
ax.text(diamond_bulk_modulus_coordinates[0], diamond_bulk_modulus_coordinates[1],
        diamond_bulk_modulus_string, color=diamond_color, horizontalalignment='center')
ax.text(beta_tin_equilibrium_volume_coordinates[0], beta_tin_equilibrium_volume_coordinates[1],
        beta_tin_equilibrium_volume_string, color='black')
ax.text(beta_tin_bulk_modulus_coordinates[0], beta_tin_bulk_modulus_coordinates[1],
        beta_tin_bulk_modulus_string, color=beta_tin_color, horizontalalignment='center')
ax.text(transition_pressure_coordinates[0], transition_pressure_coordinates[1], transition_pressure_string,
        color=tangent_line_color)

# Produce output proportional to axis labels and axes
plt.tight_layout()

# Plot output format
if write_png_file:
    plt.savefig(figure_file_name)
else:
    plt.show()

