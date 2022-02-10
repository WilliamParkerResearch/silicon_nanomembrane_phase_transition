import numpy as np
import matplotlib.pyplot as plt

# x = np.linspace(1,11,num=100)
x = np.arange(1,9)
y = 1/(x-0.25)
y2 = 1/x
plt.plot(x,y,color='blue',marker='o')
plt.plot(x,y2,color='red',marker='s')
plt.hlines(y=1,xmin=x[0],xmax=x[-1],linestyle='dashed')
plt.show()