import matplotlib as mpl

# Use TeX fonts
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.sans-serif'] = "cmr10"
mpl.rcParams['font.family'] = 'serif'

# Bands plot parameters
all_colors = 'red'
band_color = all_colors

dos_curve_color = all_colors
dos_fill_color = all_colors

dos_opacity = 0.8
fermi_level_linewidth = 0.5
fermi_level_linecolor = 'black'
fermi_level_linestyle = (0, (5, 10))
