from numpy import *
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from array import *
import linecache

dt=50
sum =0
length = 8000
data=[]

f = open("data.txt", "w")
for i in range(0,10):

	x=linecache.getline("HMC_X.dat",7000)
	xt=linecache.getline("HMC_X.dat",7000+i)
	x = x.split(" ")
	xt = xt.split(" ")
	del x[-1]
	del xt[-1]

	x = [float(a) for a in x]
	xt = [float(a) for a in xt]

	for i in range(0,len(x)):
		sum += x[i] * xt[i]
	#if(sum>0):
	data.append(sum*sqrt(5)/len(x))
	#else:
		#data.append(0)

	sum=0
for i in data:
	f.write(	str(i)	+"\n")
print(data)
plt.plot(data)
#plt.yscale('log')
plt.show()


