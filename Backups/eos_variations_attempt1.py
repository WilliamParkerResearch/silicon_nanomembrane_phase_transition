import scipy.optimize as sp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.optimize as sp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


# Parameters
number_of_volume_points = 1000
exchange_correlation = 'PBE'
chemical_formula = 'Si'
structure_names = ['$Fd\overline{3}m$', '$I4_{1}/amd$']  # Hermann-Mauguin notation for diamond & beta-Sn structures
n_atom_diamond = 8
n_atom_betasn = 8

# Conversion factors
cubic_meters_per_cubic_angstrom = 1e-30
joules_per_Rydberg = 2.1798741e-18

# Data
if exchange_correlation == 'PBE':
    figure_file_name = 'Si.PBE_1ML.EoS.png'
    volumes_sim_betasn = cubic_meters_per_cubic_angstrom*np.array([9.65158069550,10.29501940960,10.93845812370,11.58189683782,12.86877426200,14.15565168618,15.44252911440,16.72940654262])
    total_energies_strain_betasn = (joules_per_Rydberg/n_atom_betasn)*np.array([-372.66327325,-373.06805755,-373.30465661,-373.43427048,-373.50849372,-373.47259800,-373.40956222,-373.34922248])
    volumes_sim_diamond = cubic_meters_per_cubic_angstrom*np.array([14.73811721082,16.84356252521,18.94900784461,21.05445315901,23.15989847344,25.26534379280,27.37078910724])
    total_energies_strain_diamond = (joules_per_Rydberg/n_atom_betasn)*np.array([-373.29868028,-373.41432946,-373.48607830,-373.51022345,-373.48973472,-373.44785756,-373.38548468])
if exchange_correlation == 'PZ':
    figure_file_name = 'Si.PZ_1ML.EoS.png'
    volumes_sim_betasn = cubic_meters_per_cubic_angstrom*np.array([8.69553212657,9.93775100516,11.17996987587,12.42218875449,13.66440763310,14.90662650382,16.14884538240])
    total_energies_strain_betasn = (joules_per_Rydberg/n_atom_betasn)*np.array([-89.22482893,-90.30796765,-90.68395156,-90.76086191,-90.72389896,-90.65947349,-90.59691208])
    volumes_sim_diamond = cubic_meters_per_cubic_angstrom*np.array([14.08392451545,16.09591373539,18.10790294570,20.11989216564,22.13188138559,24.14387059586,26.15585981581])
    total_energies_strain_diamond = (joules_per_Rydberg/n_atom_betasn)*np.array([-90.51829730,-90.62944720,-90.70666769,-90.73095209,-90.70955754,-90.66993654,-90.60806563])
if exchange_correlation == 'SCAN':
    figure_file_name = 'Si.SCAN_1ML.EoS.png'
    volumes_sim_betasn = cubic_meters_per_cubic_angstrom*np.array([10.21619778257,11.67565460811,13.13511143364,14.59456825921,16.05402508474,17.51348191028,18.97293873584])
    total_energies_strain_betasn = (joules_per_Rydberg/n_atom_betasn)*np.array([-61.24625462,-62.21356102,-62.54424616,-62.61104424,-62.57780653,-62.51573929,-62.45101265])
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
volumes_betasn = np.linspace(volumes_sim_betasn[0], volumes_sim_betasn[-1], num=number_of_volume_points)
#
        #cell_dofree ='z' values
initial_parameters_shape_betasn = (total_energies_strain_betasn[mid(total_energies_strain_betasn)], 1e11, 3.5, volumes_sim_betasn[mid(volumes_sim_betasn)])
fit_parameters_betasn = sp.fmin(square_differences, initial_parameters_shape_betasn, args=(volumes_sim_betasn, total_energies_strain_betasn, eos), maxiter=100000)
fit_total_energies_strain_betasn = eos(fit_parameters_betasn, volumes_betasn)

# finding the negative pressure equaitons

def eosprime(p,v):              #dE/dV = -P
    dx = 1e-40
    fp = (eos(p,(v+dx)) - eos(p,v))/dx
    return fp



energyprime_diamond=eosprime(fit_parameters_diamond,volumes_diamond)
energyprime_betasn=eosprime(fit_parameters_betasn,volumes_betasn)


plt.plot(volumes_betasn,energyprime_betasn)
plt.plot(volumes_diamond,energyprime_diamond)
plt.show()
def tslope(v1,v2):
    f=(eos(fit_parameters_diamond,v2)-eos(fit_parameters_betasn,v1))/(v2-v1)
    return f


n=0
diamondidx_betaidx=[]
volume_beta_array = []
volume_diamond_array = []
slopearray = []
for i in energyprime_betasn:
    indexes = np.argwhere(np.diff(np.nan_to_num(np.sign(energyprime_diamond-i)))).flatten()
    index= np.average(indexes)
    if not np.isnan(index):
        index=int(index)

    if not np.isnan(index):
        volumes_d_b_idx = [n,index]
        # print(volumes_d_b_idx)
        volume_beta_array=np.append(volume_beta_array,volumes_betasn[volumes_d_b_idx[1]])
        volume_diamond_array=np.append(volume_diamond_array,volumes_diamond[volumes_d_b_idx[0]])

        slope=tslope(volumes_betasn[volumes_d_b_idx[1]],volumes_diamond[volumes_d_b_idx[0]])
        slopearray = np.append(slopearray,slope)
    # if index.size == 0:
    #     index=[[0]]
    # elif index.size>1:
    #     index=index[0]
    n+=1

n=0
slopediff_array=[]
for i in slopearray:
    slopediff=eosprime(fit_parameters_betasn,volume_beta_array[n])-i
    slopediff_array=np.append(slopediff_array,slopediff)
    # print(eosprime(fit_parameters_diamond,volume_diamond_array[n])-i)

    n+=1
# print(slopediff_array)

slopediff_0_idx = int(np.average(np.argwhere(np.diff(np.nan_to_num(np.sign(slopediff_array)))).flatten()))

# print(slopediff_0_idx)
tvol_diamond=volume_diamond_array[slopediff_0_idx]
tvol_betasn=volume_beta_array[slopediff_0_idx]
tpress=tslope(tvol_betasn,tvol_diamond)
tpressb=eosprime(fit_parameters_betasn,tvol_betasn)
tpressd=eosprime(fit_parameters_diamond,tvol_diamond)


print(tvol_betasn)
print(tvol_diamond)
print(tpress*1e-9)
print(eosprime(fit_parameters_diamond,tvol_diamond)*1e-9)
print(eosprime(fit_parameters_betasn,tvol_betasn)*1e-9)


vol=tvol_betasn
plt.plot(volumes_diamond,fit_total_energies_strain_diamond)
plt.plot(volumes_betasn,fit_total_energies_strain_betasn)
plt.plot(volumes_diamond,tpressb*(volumes_betasn-vol)+eos(fit_parameters_betasn,vol))
plt.show()
