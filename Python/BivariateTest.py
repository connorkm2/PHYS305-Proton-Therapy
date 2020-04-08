import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import math

def bivG(X, Y, sigma):
    Z = (1 / (2 * math.pi * sigma * sigma)) * np.exp(
        -0.5 * ((np.power(X*80, 2) / np.power(sigma, 2)) + (np.power(Y*80, 2) / np.power(sigma, 2))))
    return Z

sigma_x = 7
sigma_y = 2

x, y = np.meshgrid(np.linspace(-0.2, 0.2, 20) , np.linspace(-0.2, 0.2, 20))
#print(x)
z = bivG(x,y, sigma_x)

# for i in range(len(z)):
#     z[i] = (1 / (2 * math.pi * sigma_x * sigma_y)) * np.exp(
#         -0.5 * ((math.pow(x[i], 2) / math.pow(sigma_x, 2)) + (math.pow(y[i], 2) / math.pow(sigma_y, 2))))
#     print(z[i])

#fig = plt.figure()
ax = plt.axes(projection='3d')

ax.plot_wireframe(x, y, z*5)
ax.set_zlim(0, 0.06)

plt.show()