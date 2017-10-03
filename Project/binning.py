from numpy import *
import numpy as np
import random
from matplotlib import pyplot as plt

#data =loadtxt("HMC_Results.dat", unpack=True)
myarray = np.fromfile('HMC_Results.dat')





plt.hist(data2, bins=75)
axes = plt.gca()
#axes.set_xlim([-0.5,0.5])
plt.show()

