import matplotlib as mpl
import numpy as np

rgbcode_diamond = (60 / 255, 80 / 255, 155 / 255)
rgbcode_betasn = (156 / 255, 61 / 255, 61 / 255)
rgbcode_pressure = (61/255, 156/255, 111/255)
rgbcode_black = (44/255, 44/255, 44/255)
rgbcode_white = (235/255, 235/255, 235/255)
universal_linewidth = 1.5
axis_fontsize = 10
tick_fontsize = 12
marker_size = 50
mark_d = 's'
mark_b = 'd'
mark_p = 'o'
mark_e = '^'

cm=1/2.4
s = 9*cm
m = 14*cm
l = 19*cm


figsize_11 = np.array([1,1])
figsize_115 = np.array([1,1/1.5])
figsize_12 = np.array([1,1/2])
figsize_13 = np.array([1,1/3])
figsize_CH8 = np.array([2,1/2])
figsize_23 = np.array([1,2/3])
print(l/4)


fs_s_ch8 = tuple(s*figsize_CH8)
fs_s_11 = tuple(s*figsize_11)
fs_s_115 = tuple(s*figsize_115)
fs_s_12 = tuple(s*figsize_12)
fs_s_13 = tuple(s*figsize_13)
fs_s_23 = tuple(s*figsize_23)
fs_m_11 = tuple(m*figsize_11)
fs_m_115 = tuple(m*figsize_115)
fs_m_12 = tuple(m*figsize_12)
fs_m_13 = tuple(m*figsize_13)
fs_m_23 = tuple(m*figsize_23)
fs_l_11 = tuple(l*figsize_11)
fs_l_115 = tuple(l*figsize_115)
fs_l_12 = tuple(l*figsize_12)
fs_l_13 = tuple(l*figsize_13)
fs_l_23 = tuple(l*figsize_23)

def adjust_lightness(color, amount=0.5):
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])


rgbcode_diamond_l = adjust_lightness(rgbcode_diamond,1.3)
rgbcode_diamond_m = adjust_lightness(rgbcode_diamond,1)
rgbcode_diamond_d = adjust_lightness(rgbcode_diamond,0.7)

rgbcode_betasn_l = adjust_lightness(rgbcode_betasn,1.3)
rgbcode_betasn_m = adjust_lightness(rgbcode_betasn,1)
rgbcode_betasn_d = adjust_lightness(rgbcode_betasn,0.7)


def plot(s,name='figure.png'):
    import matplotlib.pyplot as plt
    if s==True:
        plt.show()
    else:
        plt.savefig('Figures/' + name, dpi=600, bbox_inches='tight')