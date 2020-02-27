import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

mpl.rcParams['font.serif'] = "Times"
mpl.rcParams['font.family'] = "serif"
mpl.rcParams['font.size'] = 24

exchange_correlations = ['LDA-PZ', 'GGA-PBE', 'metaGGA-SCAN']

bulk_equilibrium_volumes_per_atom = {'diamond': [19.72, 20.52, 19.61], 'beta-tin': [14.81, 15.01, 15.33]}
slab_equilibrium_volumes_per_atom = {'diamond': [20.72, 21.52, 20.61], 'beta-tin': [16.81, 17.01, 16.33]} # These values are placeholders

bulk_transition_pressures = [6.63, 8.68, 16.12]
slab_transition_pressures = [5.6, 7.0, 14.0] # These values are placeholders

horizontal_positions = np.arange(len(exchange_correlations))
bar_width = 0.35

def label_bars(rectangles):
    for rectangle in rectangles:
        height = rectangle.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rectangle.get_x() + rectangle.get_width()/2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom')


# Transition pressure plot
fig, ax = plt.subplots()
bulk_rectangles = ax.bar(horizontal_positions - bar_width/2,
                         bulk_transition_pressures, bar_width,
                         label='Bulk',
                         color='gray')
slab_rectangles = ax.bar(horizontal_positions + bar_width/2,
                         slab_transition_pressures, bar_width,
                         label='Slab',
                         color='green')
ax.set_ylabel('Transition Pressure (GPa)')
ax.set_xticks(horizontal_positions)
ax.set_xticklabels(exchange_correlations)
ax.legend(loc='upper left')

label_bars(bulk_rectangles)
label_bars(slab_rectangles)

#Plot and label experimental range
experimental_transition_pressure_minimum = 13.3
experimental_transition_pressure_maximum = 13.3
experimental_transition_pressure_mean = (experimental_transition_pressure_maximum +
                                         experimental_transition_pressure_minimum) / 2.0
experimental_transition_pressure_mean_string = '{:04.2f}'.format(experimental_transition_pressure_mean)

ax.axhline(experimental_transition_pressure_minimum, color='black',linestyle='--')
ax.axhline(experimental_transition_pressure_maximum, color='black',linestyle='--')
label_shift = 1.5*bar_width
ax.text(horizontal_positions[-1]+label_shift, experimental_transition_pressure_mean,
        'Experiment\n(Bulk)\n'+experimental_transition_pressure_mean_string,
        verticalalignment='center')

fig.tight_layout()

plt.savefig('Si.Transition_pressure.png')

# Volumes per atom diamond plot
fig, ax = plt.subplots()
bulk_rectangles = ax.bar(horizontal_positions - bar_width/2,
                         bulk_equilibrium_volumes_per_atom['diamond'], bar_width,
                         label='Bulk',
                         color='gray')
slab_rectangles = ax.bar(horizontal_positions + bar_width/2,
                         slab_equilibrium_volumes_per_atom['diamond'], bar_width,
                         label='Slab',
                         color='green')
ax.set_ylabel(r'Equilibrium Volume (Å$^3$/atom)')
ax.set_xticks(horizontal_positions)
ax.set_xticklabels(exchange_correlations)
ax.legend(loc='upper left')

label_bars(bulk_rectangles)
label_bars(slab_rectangles)

#Plot and label experimental value
diamond_experimental_lattice_constant = 5.431021 # * 10^-10 m from CODATA 2018
diamond_experimental_volume_per_atom = diamond_experimental_lattice_constant**3 / 8
diamond_experimental_volume_per_atom_string = '{:04.2f}'.format(diamond_experimental_volume_per_atom)
ax.axhline(diamond_experimental_volume_per_atom, color='black')
label_shift = 1.5*bar_width
ax.text(horizontal_positions[-1]+label_shift, diamond_experimental_volume_per_atom,
        'Experiment\n(Bulk)\n'+diamond_experimental_volume_per_atom_string,
        verticalalignment='center')

#Scale graph to be +/- 10% experimental value
factor_away_from_experiment = 0.1
ax.set_ylim([(1.0 - factor_away_from_experiment)*diamond_experimental_volume_per_atom,
             (1.0 + factor_away_from_experiment)*diamond_experimental_volume_per_atom])

fig.tight_layout()

plt.savefig('Si.Volume_per_atom.png')

