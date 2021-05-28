import scipy.optimize as sp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.optimize as sp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


# Parameters
number_of_volume_points = 100
exchange_correlation = 'PBE'
chemical_formula = 'Si'
structure_names = ['$Fd\overline{3}m$', '$I4_{1}/amd$']  # Hermann-Mauguin notation for diamond & beta-Sn structures
n_atom_diamond = 8
n_atom_betasn = 4

# Conversion factors
cubic_meters_per_cubic_angstrom = 1e-30
joules_per_Rydberg = 2.1798741e-18

# Data
if exchange_correlation == 'PBE':
    figure_file_name = 'Si.PBE_1ML.EoS.png'
    # volumes_sim_BetaSn = cubic_meters_per_cubic_angstrom*np.array([9.65158069550,10.29501940960,10.93845812370,11.58189683782,12.86877426200,14.15565168618,15.44252911440,16.72940654262])
    # total_energies_strain_BetaSn = (joules_per_Rydberg/n_atom_betasn)*np.array([-372.66327325,-373.06805755,-373.30465661,-373.43427048,-373.50849372,-373.47259800,-373.40956222,-373.34922248])
    # volumes_sim_diamond = cubic_meters_per_cubic_angstrom*np.array([14.73811721082,16.84356252521,18.94900784461,21.05445315901,23.15989847344,25.26534379280,27.37078910724])
    # total_energies_strain_diamond = (joules_per_Rydberg/n_atom_betasn)*np.array([-373.29868028,-373.41432946,-373.48607830,-373.51022345,-373.48973472,-373.44785756,-373.38548468])
    total_energies_strain_diamond = (joules_per_Rydberg / n_atom_betasn) * np.array([-373.79919010, -373.79751321, -373.79780576, -373.79253716, -373.79306892, -373.78450623, -373.78520499,-373.77376787, -373.77448566, -373.76059298, -373.76128060])
    volumes_sim_diamond = cubic_meters_per_cubic_angstrom * np.array([19.40819764559, 19.61249446291, 19.81679128023, 20.02108809756, 20.22538491488, 20.42968173220, 20.63397854952,20.83827536684, 21.04257218417, 21.24686900149, 21.45116581881])
    total_energies_strain_BetaSn = (joules_per_Rydberg / n_atom_betasn) * np.array([-186.82250887, -186.82234408, -186.82237233, -186.82186036, -186.82194978, -186.82103956, -186.82125772,-186.81986151, -186.82031125, -186.81830725, -186.81912382])
    volumes_sim_BetaSn = cubic_meters_per_cubic_angstrom * np.array([14.57700661522, 14.73044879012, 14.88389096501, 15.03733313991, 15.19077531481, 15.34421748971, 15.49765966460,15.65110183950, 15.80454401440, 15.95798618929, 16.11142836419])

if exchange_correlation == 'PZ':
    figure_file_name = 'Si.PZ_1ML.EoS.png'
    volumes_sim_BetaSn = cubic_meters_per_cubic_angstrom*np.array([8.69553212657,9.93775100516,11.17996987587,12.42218875449,13.66440763310,14.90662650382,16.14884538240])
    total_energies_strain_BetaSn = (joules_per_Rydberg/n_atom_betasn)*np.array([-89.22482893,-90.30796765,-90.68395156,-90.76086191,-90.72389896,-90.65947349,-90.59691208])
    volumes_sim_diamond = cubic_meters_per_cubic_angstrom*np.array([14.08392451545,16.09591373539,18.10790294570,20.11989216564,22.13188138559,24.14387059586,26.15585981581])
    total_energies_strain_diamond = (joules_per_Rydberg/n_atom_betasn)*np.array([-90.51829730,-90.62944720,-90.70666769,-90.73095209,-90.70955754,-90.66993654,-90.60806563])
if exchange_correlation == 'SCAN':
    figure_file_name = 'Si.SCAN_1ML.EoS.png'
    volumes_sim_BetaSn = cubic_meters_per_cubic_angstrom*np.array([10.21619778257,11.67565460811,13.13511143364,14.59456825921,16.05402508474,17.51348191028,18.97293873584])
    total_energies_strain_BetaSn = (joules_per_Rydberg/n_atom_betasn)*np.array([-61.24625462,-62.21356102,-62.54424616,-62.61104424,-62.57780653,-62.51573929,-62.45101265])
    volumes_sim_diamond = cubic_meters_per_cubic_angstrom*np.array([13.94338719066,15.93529964919,17.92721209824,19.91912455677,21.91103701526,23.90294946430,25.89486192283])
    total_energies_strain_diamond = (joules_per_Rydberg/n_atom_diamond)*np.array([-62.34036437,-62.49221543,-62.53738742,-62.53699016,-62.50937811,-62.43660083,-62.35068451])
def mid(x):
    return int((len(x)-1)/2)


def square_differences(p, x, y, f):
    return np.power(f(p, x) - y, 2).sum()


def eos(p, v):        #eos
    kk = p[2]-1.0
    return (p[0] + (p[1]*p[3]*(((1.0/(p[2]*kk))*np.power((v/p[3]), (-kk)))+(v/(p[2]*p[3]))-(1.0/kk))))


# eos Equation of State Information
#           parameters = (E0, K0, K0', V0)

#   Diamond structure calculations
initial_parameters_diamond = (total_energies_strain_diamond[mid(total_energies_strain_diamond)], 1e11, 3.5, volumes_sim_diamond[mid(volumes_sim_diamond)])
fit_parameters_diamond = sp.fmin(square_differences, initial_parameters_diamond, args=(volumes_sim_diamond, total_energies_strain_diamond, eos), maxiter=100000)
volumes_diamond = np.linspace(volumes_sim_diamond[0], volumes_sim_diamond[-1], num=number_of_volume_points)
fit_total_energies_strain_diamond = eos(fit_parameters_diamond, volumes_diamond)

#     Beta-Sn structure calculations
volumes_BetaSn = np.linspace(volumes_sim_BetaSn[0], volumes_sim_BetaSn[-1], num=number_of_volume_points)
#
        #cell_dofree ='z' values
initial_parameters_shape_BetaSn = (total_energies_strain_BetaSn[mid(total_energies_strain_BetaSn)], 1e11, 3.5, volumes_sim_BetaSn[mid(volumes_sim_BetaSn)])
fit_parameters_shape_BetaSn = sp.fmin(square_differences, initial_parameters_shape_BetaSn, args=(volumes_sim_BetaSn, total_energies_strain_BetaSn, eos), maxiter=100000)
fit_total_energies_strain_BetaSn = eos(fit_parameters_shape_BetaSn, volumes_BetaSn)

# finding the negative pressure equaitons

def eosprime(p,v):              #dE/dV = -P
    dx = 1e-40
    fp = (eos(p,(v+dx)) - eos(p,v))/dx
    return fp



# energyprime_diamond=eosprime(fit_parameters_diamond,volumes_diamond)
# energyprime_beta=eosprime(fit_parameters_shape_BetaSn,volumes_BetaSn)


# plt.plot(volumes_BetaSn,energyprime_beta)
# plt.plot(volumes_diamond,energyprime_diamond)
# plt.show()

def tslope(v1,v2):
    f=(eos(fit_parameters_diamond,v2)-eos(fit_parameters_shape_BetaSn,v1))/(v2-v1)
    return f

#finds volumes where the tangent line of the beta eos goes from 2 intersection points with diamond eos to 0

def tangent_beta(vol0,vol):
    t = eosprime(fit_parameters_shape_BetaSn,vol0)*(vol-vol0)+eos(fit_parameters_shape_BetaSn,vol0)
    return t


def intersection_idx(equation):     #prints lowerbound index in wihch the value is found within it and the index above it
    i = np.argwhere(np.diff(np.nan_to_num(np.sign(equation)))).flatten()
    return i

pos_tvol_diamond1,pos_tvol_diamond2 = np.amin(volumes_diamond), np.amax(volumes_diamond)
pos_tvol_beta1,pos_tvol_beta2 = np.amin(volumes_BetaSn), np.amax(volumes_BetaSn)
slopediff=1e10
slope_precision = 0.01*1e9
numvol=1000
while slopediff >= slope_precision:
    volumes_diamond_loop = np.linspace(pos_tvol_diamond1,pos_tvol_diamond2,numvol)
    volumes_BetaSn_loop = np.linspace(pos_tvol_beta1,pos_tvol_beta2,numvol)
    intersection_num_array=[]
    for i in volumes_BetaSn_loop:
        intersection_equations=tangent_beta(i,volumes_diamond_loop)-eos(fit_parameters_diamond,volumes_diamond_loop)  #used to solve: tangentlines - E_diamond =0
        intersection_points= intersection_idx(intersection_equations) #finds intersection point
        num_int = len(intersection_points)
        intersection_num_array=np.append(intersection_num_array,num_int)
    # print(np.diff(intersection_num_array))

    turnpoint_idx_beta_1, turnpoint_idx_beta_2 = int(np.where(np.diff(intersection_num_array)==2)[0][0])-10, int(np.where(np.diff(intersection_num_array)==2)[0][0])+10
    print(turnpoint_idx_beta_1,turnpoint_idx_beta_2)
    pos_tvol_beta1, pos_tvol_beta2 = volumes_BetaSn_loop[turnpoint_idx_beta_1],volumes_BetaSn_loop[turnpoint_idx_beta_2]

    intersection_equations1 = tangent_beta(pos_tvol_beta1, volumes_diamond_loop) - eos(fit_parameters_diamond, volumes_diamond_loop)
    intersection_equations2 = tangent_beta(pos_tvol_beta2, volumes_diamond_loop) - eos(fit_parameters_diamond, volumes_diamond_loop)
    intersection_points_idx_diamond_1 = intersection_idx(intersection_equations1)
    intersection_points_idx_diamond_2 = intersection_idx(intersection_equations2)
    print(intersection_points_idx_diamond_1)
    print(intersection_points_idx_diamond_2)

    if int(len(intersection_points_idx_diamond_1))==2:
        pos_tvol_diamond1 = volumes_diamond_loop[intersection_points_idx_diamond_1[0]]
        pos_tvol_diamond2 = volumes_diamond_loop[intersection_points_idx_diamond_1[1]+1]
    if int(len(intersection_points_idx_diamond_2))==2:
        pos_tvol_diamond1 = volumes_diamond_loop[intersection_points_idx_diamond_2[0]]
        pos_tvol_diamond2 = volumes_diamond_loop[intersection_points_idx_diamond_2[1]+1]


    # tvol_diamond1= volumes_diamond_loop[intersection_points_idx_diamond_1[0]+10]
    # tvol_diamond2= volumes_diamond_loop[intersection_points_idx_diamond_1[0]-10]


    slope1, slope2 = eosprime(fit_parameters_shape_BetaSn,pos_tvol_beta1),eosprime(fit_parameters_shape_BetaSn,pos_tvol_beta2)

    slopediff = np.absolute(slope2-slope1)
    print(pos_tvol_beta1,pos_tvol_beta2)
    print(pos_tvol_beta1-pos_tvol_beta2)
    print(pos_tvol_diamond1,pos_tvol_diamond2)
    print(pos_tvol_diamond1-pos_tvol_diamond2)
    print(eosprime(fit_parameters_shape_BetaSn,pos_tvol_beta1)*1e-9,eosprime(fit_parameters_shape_BetaSn,pos_tvol_beta2)*1e-9)
    print(slopediff*1e-9)

    # plt.plot(volumes_diamond_loop,eos(fit_parameters_shape_BetaSn,volumes_diamond_loop))
    # plt.plot(volumes_diamond_loop,tangent_beta(pos_tvol_beta1,volumes_diamond_loop))
    # plt.plot(volumes_diamond_loop,tangent_beta(pos_tvol_beta2,volumes_diamond_loop))
    # plt.show()
volumes_range = np.append(volumes_BetaSn,volumes_diamond)
plt.plot(volumes_BetaSn,fit_total_energies_strain_BetaSn)
plt.plot(volumes_diamond,fit_total_energies_strain_diamond)
plt.plot(volumes_range,tangent_beta(pos_tvol_beta1,volumes_range))
plt.plot(volumes_range,tangent_beta(pos_tvol_beta2,volumes_range))
plt.show()
