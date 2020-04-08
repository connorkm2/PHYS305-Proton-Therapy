
## README ##
#
# This is a simple script for quickly plotting all useful figures from the output csv files.
# The figures generated here are to be properly formatted before use in the final report.
#
# To use this script update the paths for the csv files and change the slice number if you are using
# a different slice number in the main simulation code.
#
# Will need to add more figures, for example the individual Bragg peaks plot.
#

import os
import matplotlib.pyplot as plt
import numpy as np
os.chdir("../")

data = np.genfromtxt("15_zSlice.csv", delimiter=",")

plt.figure(1)
plt.imshow(data, cmap='hot', interpolation='nearest')
plt.gray()
plt.xlabel("x")
plt.ylabel("y")
plt.colorbar()


# Dose depth distribution
data = np.genfromtxt("energy_hist.csv", delimiter=",")

plt.figure(2)
plt.plot(data[:, 4], data[:, 5])
plt.xlim(0.22, 0.52)
plt.xlabel("Depth [m]")
plt.ylabel("Energy")

# Dose profile
data = np.genfromtxt("dose_profile.csv", delimiter=",")
plt.figure(3)
plt.title("Dose Profile")
plt.plot(data[:, 0], data[:, 1])
plt.ylabel("Energy")

# Detector
data = np.genfromtxt("det_theta_zx.csv", delimiter=",")
plt.figure(4)
plt.title("Detector angles")
plt.plot(data[:, 1], data[:, 2])
plt.xlabel("Theta")

plt.show()
