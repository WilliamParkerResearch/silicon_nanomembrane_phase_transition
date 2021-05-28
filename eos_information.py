import scipy.optimize as sp
import numpy as np
import matplotlib.pyplot as plt
from ReferenceFiles.FunctionDefinitions import *
from ReferenceFiles.unit_conversions import *

# Parameters
number_of_volume_points = 10000
pressure_extreme = 50*11.7e10
exchange_correlation = 'PBE'
N_ML = '2'

directoryofdata='DataFolder'+'.'+exchange_correlation+'.'+'Data_'+N_ML+'L'
exec(f'from {directoryofdata} import *')


def eos(p, v):        #eos
    kk = p[2]-1.0
    return (p[0] + (p[1]*p[3]*(((1.0/(p[2]*kk))*np.power((v/p[3]), (-kk)))+(v/(p[2]*p[3]))-(1.0/kk))))


def eosprime(p,v):              #dE/dV = -P
    dx = 1e-31
    fp = (eos(p,(v+dx)) - eos(p,v))/dx
    return fp


def tslope(v1,v2):
    f=(eos(fit_parameters_diamond,v2)-eos(fit_parameters_betasn,v1))/(v2-v1)
    return f

#finds volumes where the tangent line of the beta eos goes from 2 intersection points with diamond eos to 0

def tangent_beta(vol0,vol):
    t = eosprime(fit_parameters_betasn,vol0)*(vol-vol0)+eos(fit_parameters_betasn,vol0)
    return t


def intersection_idx(equation):     #prints lowerbound index in wihch the value is found within it and the index above it
    i = np.argwhere(np.diff(np.nan_to_num(np.sign(equation)))).flatten()
    return i


# eos Equation of State Information
#           parameters = (E0, K0, K0', V0)

#   Diamond structure calculations
initial_parameters_diamond = (total_energies_strain_diamond[mid(total_energies_strain_diamond)], 1e11, 3.5, volumes_sim_diamond[mid(volumes_sim_diamond)])
fit_parameters_diamond = sp.fmin(square_differences, initial_parameters_diamond, args=(volumes_sim_diamond, total_energies_strain_diamond, eos), maxiter=100000)
volumes_diamond = np.linspace(volumes_sim_diamond[0], volumes_sim_diamond[-1], num=number_of_volume_points)
fit_total_energies_strain_diamond = eos(fit_parameters_diamond, volumes_diamond)

#     Beta-Sn structure calculations
initial_parameters_shape_betasn = (total_energies_strain_betasn[mid(total_energies_strain_betasn)], 1e11, 3.5, volumes_sim_betasn[mid(volumes_sim_betasn)])
fit_parameters_betasn = sp.fmin(square_differences, initial_parameters_shape_betasn, args=(volumes_sim_betasn, total_energies_strain_betasn, eos), maxiter=100000)
volumes_betasn = np.linspace(volumes_sim_betasn[0], volumes_sim_betasn[-1], num=number_of_volume_points)
fit_total_energies_strain_betasn = eos(fit_parameters_betasn, volumes_betasn)


pos_tvol_diamond1,pos_tvol_diamond2 = np.amin(volumes_diamond), np.amax(volumes_diamond)
pos_tvol_beta1,pos_tvol_beta2 = np.amin(volumes_betasn), np.amax(volumes_betasn)
slopediff=1e10
slope_precision = 1e-3
numvol=1000
while slopediff >= slope_precision:
    volumes_diamond_loop = np.linspace(pos_tvol_diamond1,pos_tvol_diamond2,numvol)
    volumes_betasn_loop = np.linspace(pos_tvol_beta1,pos_tvol_beta2,numvol)
    intersection_num_array=[]
    for i in volumes_betasn_loop:
        intersection_equations=tangent_beta(i,volumes_diamond_loop)-eos(fit_parameters_diamond,volumes_diamond_loop)  #used to solve: tangentlines - E_diamond =0
        intersection_points = intersection_idx(intersection_equations) #finds intersection point
        num_int = len(intersection_points)
        intersection_num_array=np.append(intersection_num_array,num_int)


    turnpoint_idx_beta_1, turnpoint_idx_beta_2 = int(np.where(np.diff(intersection_num_array)==2)[0][0])-10, int(np.where(np.diff(intersection_num_array)==2)[0][0])+10
    pos_tvol_beta1, pos_tvol_beta2 = volumes_betasn_loop[turnpoint_idx_beta_1],volumes_betasn_loop[turnpoint_idx_beta_2]

    intersection_equations1 = tangent_beta(pos_tvol_beta1, volumes_diamond_loop) - eos(fit_parameters_diamond, volumes_diamond_loop)
    intersection_equations2 = tangent_beta(pos_tvol_beta2, volumes_diamond_loop) - eos(fit_parameters_diamond, volumes_diamond_loop)
    intersection_points_idx_diamond_1 = intersection_idx(intersection_equations1)
    intersection_points_idx_diamond_2 = intersection_idx(intersection_equations2)


    if int(len(intersection_points_idx_diamond_1))==2:
        pos_tvol_diamond1 = volumes_diamond_loop[intersection_points_idx_diamond_1[0]]
        pos_tvol_diamond2 = volumes_diamond_loop[intersection_points_idx_diamond_1[1]+1]
    if int(len(intersection_points_idx_diamond_2))==2:
        pos_tvol_diamond1 = volumes_diamond_loop[intersection_points_idx_diamond_2[0]]
        pos_tvol_diamond2 = volumes_diamond_loop[intersection_points_idx_diamond_2[1]+1]


    slope1, slope2 = eosprime(fit_parameters_betasn,pos_tvol_beta1),eosprime(fit_parameters_betasn,pos_tvol_beta2)

    slopediff = np.absolute(slope2-slope1)
print('volt_b:',pos_tvol_beta1*1e30,'-',pos_tvol_beta2*1e30)
print('delta_volt_b:',(pos_tvol_beta1-pos_tvol_beta2)*1e30)
print('volt_d:',pos_tvol_diamond1*1e30,'-',pos_tvol_diamond2*1e30)
print('delta volt_b:',(pos_tvol_diamond1-pos_tvol_diamond2)*1e30)
print('pt:',eosprime(fit_parameters_betasn,pos_tvol_beta1)*1e-9,'-',eosprime(fit_parameters_betasn,pos_tvol_beta2)*1e-9)
print('delta pt',slopediff*1e-9)


volumes_range = np.append(volumes_betasn,volumes_diamond)
plt.plot(volumes_betasn,fit_total_energies_strain_betasn)
plt.plot(volumes_diamond,fit_total_energies_strain_diamond)
plt.plot(volumes_range,tangent_beta(pos_tvol_beta1,volumes_range))
plt.plot(volumes_range,tangent_beta(pos_tvol_beta2,volumes_range))
plt.show()
