import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import make_interp_spline
from ReferenceFiles.format_charts import *

# mpl.rcParams['font.serif'] = "Times"
# mpl.rcParams['font.family'] = "Georgia"
mpl.rcParams['font.size'] = 20
b_pt_color='limegreen'
s_pt_color='darkseagreen'
ss_pt_color='green'

s_layers=np.array([1,2,3])
s_pt=np.array([-1.35,3.8,4.86])

##PBE
layers = np.array([1,2,3,4,5,6])
bulkpt = 8.44
pt = np.array([-1.46,0.06,2.63,4.27,4.97,6.14])
cellad=np.array([10.32618687240,10.22701581000,10.27866740500,10.28030461780000,10.28030461780000])
cellcd=np.array([7.28118500750800681633,8.12285311560294625329,8.53192071244818663869,8.74451990906822779714,8.87689169493783859552])
cellab = np.array([9.01462787842,9.10522715358,9.05992751600,9.05992751600,9.05992751600])
cellcb = np.array([5.99492789376121182371,7.14686575725938618693,7.44464338107961951400,7.61418761059102109311,7.69251098294275038470])
cellad_bulk = 10.33190348
cellcd_bulk = 10.33190348
cellab_bulk = 9.060447249
cellcb_bulk = 5.043749108


# # # #SCAN
# layers = np.array([1,2,3,4])
# bulkpt = 16.07
# pt = np.array([0.29,4.04,8.72,9.89])
# cellad_bulk = 10.33196444
# cellcd_bulk = 10.33196444
# cellab_bulk = 9.059927516
# cellcb_bulk = 5.046060010
# cellad=np.array([10.3819705950000000,10.2786674050000000,10.22727406797])
# cellcd=np.array([7.13797304985846315391,7.87584247489052551500,8.28044250935143524922])
# cellab = np.array([9.0146278784200,9.05992751600,9.0146278784200])
# cellcb = np.array([5.79624949631693115740,7.03454880545041841523,7.48240357093339795504])







layers_new = np.linspace(layers.min(),layers.max(),300)
s_layers_new = np.linspace(s_layers.min(),s_layers.max(),300)

spl = make_interp_spline(layers,pt,k=2)
pt_smooth = spl(layers_new)

sspl = make_interp_spline(s_layers,s_pt,k=2)
s_pt_smooth = sspl(s_layers_new)

plt.plot(layers_new,pt_smooth,label='Original Relaxed',color=s_pt_color,linestyle='dotted')
plt.scatter(layers,pt,color=s_pt_color)
plt.plot(s_layers_new,s_pt_smooth,label='More Stable',color=ss_pt_color,linestyle='dotted')
plt.scatter(s_layers,s_pt,color=ss_pt_color)
# plt.title('Transition Pressure as Layers of Membranes are Added for metaGGA-SCAN')
plt.xlabel(r'$N_{\rm ML}$')
plt.xticks(np.append(layers,np.array([layers[-1]+1])))
plt.ylabel(r"$P_{\rm T}$ (GPa)")
plt.ylim(-2.5,9)
plt.axhline(y=bulkpt,color=b_pt_color,linestyle='--',label='Bulk Value Pt')
plt.axhline(y=0,color='grey')
plt.text(1,bulkpt-0.5,'Bulk',color=b_pt_color)
plt.text(1,pt[0]-0.75,'Nanomembrane',color=ss_pt_color)
plt.text(s_layers[-2]-0.25,s_pt[-1]+0.50,'More Stable Structure',color=ss_pt_color,size=15)
plt.text(s_layers[-1]+0.25,pt[2]-0.75,'Original Relaxed Structure',color=s_pt_color,size=15)


plt.show()

# plt.plot(layers,cellcb,marker='o', color = "r",label='cell_c_betatin',linestyle='dotted')
# plt.axhline(y=cellcb_bulk,color='r',linestyle='--',label='Bulk Value cellcb')
# plt.text(layers[-1],cellcb[-1]-0.5,'Beta-Tin',color='r')
# plt.plot(layers,cellcd,marker='o', color = "b",label='cell_c_diamond',linestyle='dotted')
# plt.axhline(y=cellcd_bulk,color='b',linestyle='--',label='Bulk Value cellcb')
# plt.text(layers[-1],cellcd[-1]-0.5,'Diamond',color='b')
# plt.title('Out of Plane Parameter as Layers of Membranes are Added for metaGGA-SCAN')
# plt.xlabel("Cubic Monolayer Number")
# plt.xticks(np.append(layers,np.array([layers[-1]+1])))
# plt.ylabel("Out of PLane Lattice Parameter (bohr)")
# plt.show()
#
#
#
# plt.plot(layers,cellab,marker='o',color="r",label='cell_a_betatin',linestyle='dotted')
# plt.axhline(y=cellab_bulk,color='r',linestyle='--',label='Bulk Value cellcb')
# plt.text(layers[-1],cellab[-1]-0.8*(cellab.max()-cellab.min()),'Beta-Tin',color='r')
# plt.plot(layers,cellad,marker='o',color="b",label='cell_a_diamond',linestyle='dotted')
# plt.axhline(y=cellad_bulk,color='b',linestyle='--',label='Bulk Value cellcb')
# plt.text(layers[-1],cellad[-1]-0.8*(cellad.max()-cellad.min()),'Diamond',color='b')
# plt.title('In Plane Lattice Parameter as Layers of Membranes are Added for metaGGA-SCAN')
# plt.xlabel("Cubic Monolayer Number")
# plt.xticks(np.append(layers,np.array([layers[-1]+1])))
# plt.ylabel("In Plane Lattice Parameter (bohr)")
# plt.show()