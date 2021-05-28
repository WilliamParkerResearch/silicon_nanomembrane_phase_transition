import scipy.optimize as sp
import numpy as np
import matplotlib.pyplot as plt
from ReferenceFiles.FunctionDefinitions import *
from ReferenceFiles.unit_conversions import *

# Parameters
number_of_volume_points = 10000
pressure_extreme = 50*11.7e10
exchange_correlation = 'PBE'
N_ML = '0'

directoryofdata='DataFolder'+'.'+exchange_correlation+'.'+'eos'+'.'+'Data_'+N_ML+'L'
exec(f'from {directoryofdata} import *')


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
initial_parameters_shape_betasn = (total_energies_strain_betasn[mid(total_energies_strain_betasn)], 1e11, 3.5, volumes_sim_betasn[mid(volumes_sim_betasn)])
fit_parameters_betasn = sp.fmin(square_differences, initial_parameters_shape_betasn, args=(volumes_sim_betasn, total_energies_strain_betasn, eos), maxiter=100000)
volumes_betasn = np.linspace(volumes_sim_betasn[0], volumes_sim_betasn[-1], num=number_of_volume_points)
fit_total_energies_strain_betasn = eos(fit_parameters_betasn, volumes_betasn)

numvol = 100
vol_d_del_0 = 1e-2
vol_b_del_0 = 1e-3
pt_del_0 = 1e-3
vol_d_del = 1
vol_b_del = 1
pt_del = 1
vol_d1,vol_d2 = np.amin(volumes_diamond), np.amax(volumes_diamond)
vol_b1,vol_b2 = np.amin(volumes_betasn), np.amax(volumes_betasn)

volumes_d_loop = np.linspace(vol_b1 + (vol_b2 - vol_b1) / 3, vol_d2 + (vol_b2 - vol_b1), numvol)
volumes_b_loop = np.linspace(vol_b1, vol_b2, numvol)

# while vol_d_del >= vol_d_del_0:
while vol_d_del >= vol_d_del_0:
    while pt_del >= pt_del_0:
        eos_d_loop =eos(fit_parameters_diamond,volumes_d_loop)
        eos_d_loop_prime = np.gradient(eos_d_loop,volumes_b_loop)
        eos_b_loop = eos(fit_parameters_betasn,volumes_b_loop)
        eos_b_loop_prime = np.gradient(eos_b_loop,volumes_b_loop)


        def tangent_beta_loop(volb_idx,vol):
            t = eos_b_loop_prime[volb_idx]*(vol-volumes_b_loop[volb_idx])+eos(fit_parameters_betasn,volumes_b_loop[volb_idx])
            return t


        #finds volumes where the tangent line of the beta eos goes from 2 intersection points with diamond eos to 0
        def intersection_idx(equation):     #prints lowerbound index in wihch the value is found within it and the index above it
            int = np.diff(np.nan_to_num(np.sign(equation)))
            return int


        for i in range(numvol):
            plt.plot(volumes_d_loop, tangent_beta_loop(i,volumes_d_loop),color='cyan')
            intersection_equations = tangent_beta_loop(i,volumes_d_loop)-eos(fit_parameters_diamond,volumes_d_loop)
            intersection_points = intersection_idx(intersection_equations)  # finds intersection point
            print(intersection_points)
            if not np.all((intersection_points == 0)):
                vol_d_min_idx = np.where(intersection_points==2)[0][0]
                vol_d_max_idx = np.where(intersection_points==-2)[0][0]+1
                vol_d1 = volumes_d_loop[vol_d_min_idx]
                vol_d2 = volumes_d_loop[vol_d_max_idx]
                vol_d_del = np.abs(vol_d2-vol_d1)*1e30
                vol_b1 = volumes_b_loop[i-1]
                vol_b2 = volumes_b_loop[i]
                vol_b_del = np.abs(vol_b2-vol_b1)*1e30
                pt1 = eos_b_loop_prime[i-1]*1e-9
                pt2 = eos_b_loop_prime[i]*1e-9
                pt_del = np.abs(pt2-pt1)
                volumes_d_loop = np.linspace(vol_d1, vol_d2, numvol)
                volumes_b_loop = np.linspace(vol_b1, vol_b2, numvol)
                break
        print(vol_d1,'-',vol_d2,':',vol_d_del)
        print(vol_b1,'-',vol_b2,':',vol_b_del)
        print(pt1,'-',pt2,':',pt_del)

plt.plot(volumes_b_loop,eos(fit_parameters_betasn,volumes_b_loop),color='black')
plt.plot(volumes_d_loop,eos(fit_parameters_diamond,volumes_d_loop),color='black')
plt.show()