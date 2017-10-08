from numpy import *
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes

#data =loadtxt("HMC_Results.dat", unpack=True)
myarray = np.loadtxt('HMC_Results.dat',delimiter = ',',unpack = True)

#print(myarray)


out = plt.hist(myarray, bins=50,normed =1)
axes = plt.gca()
plt.xlabel('x')
plt.ylabel('|ψ|^2')
plt.title('Wave Function for a Harmonic Oscillator')
print(sum(out[0][1:49]*diff(out[1][1:50])))
plt.show()

