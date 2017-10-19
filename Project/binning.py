from numpy import *
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
import statistics 
from array import *

#data =loadtxt("HMC_Results.dat", unpack=True)
myarray = np.loadtxt('HMC_Results.dat',delimiter = ',',unpack = True)


#print(myarray)
nbins=50


#y,binEdges = np.histogram(myarray, bins=nbins)
#axes = plt.gca()

#bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
#menStd     = np.sqrt(y)
#width      = 0.05
#plt.bar(bincenters, y, width=width, color='r', yerr=menStd)
print(sum(myarray)/len(myarray))

x = np.linspace(-2,2,100) # 100 linearly spaced numbers
y = (1/(3.141**0.5))*np.exp(-x**2)
out1 = plt.plot(x,y)
out = plt.hist(myarray, bins=50,normed=1)
axes = plt.gca()
#axes.set_xlim(2,-2)

plt.xlabel('x')
plt.ylabel('|ψ|^2')
plt.title('Wave Function for a Harmonic Oscillator')
#print(sum(y[0][1:49]*diff(y[1][1:50])))
plt.show()

