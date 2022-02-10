import numpy as np
import matplotlib.pyplot as plt
from elastic_information import elastic_properties
from ReferenceFiles.plot_formating import *


max_nml = 8

showplot=True
adj_l = 1.2
# n_layers = np.arange(max_nml+1)
n_layers = np.array([0,1,2,3,4,5,7,8])
poisson_diamond = np.array([])
youngs_diamond = np.array([])
poisson_betasn = np.array([])
youngs_betasn = np.array([])

for i in n_layers:
    print(i)
    if i == 0:
        pd0,pb0,yd0,yb0 = elastic_properties(i,folder='elastic')
    else:
        pd,pb,yd,yb = elastic_properties(i,folder='elastic')
        poisson_diamond = np.append(poisson_diamond,pd)
        youngs_diamond = np.append(youngs_diamond,yd)
        poisson_betasn = np.append(poisson_betasn,pb)
        youngs_betasn = np.append(youngs_betasn,yb)



# Poissons Ratio
fig = plt.figure(figsize=fs_m_12)
plt.plot(n_layers[1:],poisson_diamond,marker=mark_d,color=rgbcode_diamond,linewidth=universal_linewidth)
plt.hlines(y=pd0,xmin=n_layers[1],xmax=max_nml,color=adjust_lightness(rgbcode_diamond,adj_l),linewidth=universal_linewidth,linestyle='dashed')

plt.plot(n_layers[1:],poisson_betasn,marker=mark_b,color=rgbcode_betasn,linewidth=universal_linewidth)
plt.hlines(y=pb0,xmin=n_layers[1],xmax=max_nml,color=adjust_lightness(rgbcode_betasn,adj_l),linewidth=universal_linewidth,linestyle='dashed')

plt.xlabel('N-ML',fontsize=axis_fontsize)
plt.ylabel(r'$\nu$',fontsize=axis_fontsize)
plt.xticks(n_layers[1:],fontsize=tick_fontsize)
plt.yticks(fontsize=tick_fontsize)
plot(showplot,'nml_poissons_ratio')


# Youngs Modulus
fig = plt.figure(figsize=fs_m_12)
plt.plot(n_layers[1:],youngs_diamond*1e-9,marker=mark_d,color=rgbcode_diamond,linewidth=universal_linewidth)
plt.hlines(y=yd0*1e-9,xmin=n_layers[1],xmax=max_nml,color=adjust_lightness(rgbcode_diamond,adj_l),linewidth=universal_linewidth,linestyle='dashed')

plt.plot(n_layers[1:],youngs_betasn*1e-9,marker=mark_b,color=rgbcode_betasn,linewidth=universal_linewidth)
plt.hlines(y=yb0*1e-9,xmin=n_layers[1],xmax=max_nml,color=adjust_lightness(rgbcode_betasn,adj_l),linewidth=universal_linewidth,linestyle='dashed')

plt.xlabel('N-ML',fontsize=axis_fontsize)
plt.ylabel('E (GPa)',fontsize=axis_fontsize)
plt.xticks(n_layers[1:],fontsize=tick_fontsize)
plt.yticks(fontsize=tick_fontsize)
plot(showplot,'nml_youngs_modulus')

