import numpy as np


# Function Definitions

def round(n, decimals=0):
    multiplier = 10 ** decimals
    return int(n * multiplier) / multiplier


def murnaghan(p, v):
    kk = p[2]-1.0
    return (p[0] + (p[1]*p[3]*(((1.0/(p[2]*kk))*np.power((v/p[3]), (-kk)))+(v/(p[2]*p[3]))-(1.0/kk))))


def square_differences(p, x, y, f):
    return np.power(f(p, x) - y, 2).sum()


def mid(x):
    return int((len(x)-1)/2)


def matchmaker(a, b):
    list_a = list(a)
    list_b = list(b)

    matches = np.intersect1d(a, b)
    index_a = []
    for x in matches:
        index_a.append(list_a.index(x))
    index_b = []
    for x in matches:
        index_b.append(list_b.index(x))
    return [index_a, index_b]


def enthal_murn(a, p):
    e0 = a[0]
    k0 = a[1]
    k0prime = a[2]
    v0 = a[3]
    if (1+p*(k0prime/k0)).any() > 0 and k0 != 0 and k0prime != 0:
        v = v0*np.power((1 + p * (k0prime / k0)), (-1.0 / k0prime))
    enthalpy = e0
    enthalpy += (k0 * v0 /(k0prime * (k0prime - 1.0))) * np.power(v/v0, (1 - k0prime))
    enthalpy += (k0 * v0 / k0prime) * (v/v0)
    enthalpy += -(k0 * v0 /(k0prime - 1))
    enthalpy += p*v
    return enthalpy


def transition_volume(a, p):
    e0 = a[0]
    k0 = a[1]
    k0prime = a[2]
    v0 = a[3]
    v = v0 * np.power((1 + p * (k0prime / k0)), (-1.0 / k0prime))
    return v


# Finds nearest value of array input to the value input
def find_nearest(array, value):
    array = np.asarray(array)
    idx = np.argmin((np.abs(array - value)))
    return [array[idx], idx]