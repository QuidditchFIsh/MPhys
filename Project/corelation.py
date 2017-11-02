from numpy import *
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from array import *
import linecache

sum =0;sum1=0;sum2=0
data=[]
data1=[]

f = open("data.txt", "w")
g = open("logData.txt","w")
for j in range(0,9000):
	x  =linecache.getline("HMC_X.dat",90000+j)
	xt =linecache.getline("HMC_X.dat",90000+1+j)
	xt1=linecache.getline("HMC_X.dat",90000+2+j)

	x   = x.split(" ")
	xt  = xt.split(" ")
	xt1 = xt1.split(" ")

	del x  [-1]
	del xt [-1]
	del xt1[-1]

	x   = [float(a) for a in x  ]
	xt  = [float(b) for b in xt ]
	xt1 = [float(c) for c in xt1]
	
	for i in range(0,len(x)):
		sum  += x[i] * xt1[i]
		sum1 += x[i] * xt [i]
		#print(str(x[i]) + " " + str(xt[i]) + " " + str(xt1[i]))

	#if(sum>0):
#	sum  = sum/len(x)
#	sum1 = sum1/len(x)
	sum2 += log(sum1/sum)
	if(j>1000):
		data.append(sum2/(j))
#else:
#data.append(0)
	#print(log(sum1/sum))

print(sum2/9000)
#plt.yscale('log')
plt.plot(data)
#plt.yscale('log')
plt.show()


