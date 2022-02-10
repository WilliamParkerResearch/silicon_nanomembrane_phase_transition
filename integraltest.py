import numpy as np
import matplotlib.pyplot as plt

x = np.arange(-3,3)
y = x

print(x)


def integral(x,y):
    delta_x = np.diff(x)
    zero_idx = np.argsort(np.abs(x))[0]
    y_0 = y[zero_idx]
    x_0 = x[zero_idx]

    area_array = np.array([])
    for i in np.arange(len(delta_x)):
        idx_bounds = np.sort(np.array([i,zero_idx]))
        idxs = np.arange(idx_bounds[0],idx_bounds[-1]+1,1)
        if i < zero_idx:
            area = -np.sum(y[idxs]*delta_x[idxs])
        else:
            area = np.sum(y[idxs]*delta_x[idxs])
        area_array = np.append(area_array,area)

    return area_array

integral(x,y)

plt.plot(x,y)
plt.plot(x[:-1],integral(x,y))
plt.show()