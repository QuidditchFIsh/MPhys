from numpy import *
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from array import *
import linecache
import statistics 

sum =0;sum1=0;sum2=0;sum3=0;
data=[]
data1=[]
data2=[]
data3=[]
data4=[]

f = open("data.txt", "w")
g = open("logData.txt","w")

t_range=50
data_range=10
for k in range(0,t_range):

	data1.append([])
	data3.append([])

	for j in range(0,data_range):
		x  =linecache.getline("HMC_X.dat",10000)
		xt =linecache.getline("HMC_X.dat",10000+j+k+1)
		xt1 =linecache.getline("HMC_X.dat",10000+j+k)

		x   = x.split(" ")
		xt  = xt.split(" ")
		xt1 = xt1.split(" ")

		del x   [-1]
		del xt  [-1]
		del xt1 [-1]

		x    = [float(a) for a in x   ]
		xt   = [float(b) for b in xt  ]
		xt1  = [float(c) for c in xt1 ]
	
		for i in range(0,len(x)):
			sum  += x[i] * xt1[i]
			sum1 += x[i] * xt [i]

		sum3 += sum/len(x)
		sum2 += sum1/len(x)

		data1[k].append(sum)
		data3[k].append(sum1)

		sum1=0;sum=0


	data.append((sum2/sum3))

	sum1=0
	sum2=0

for i in range(0,t_range):

	data2.append(np.std(data1[i]))
	data4.append(np.std(data3[i]))

#print(data)
#x = np.linspace(0, t_range,num = t_range)
#plt.errorbar(x,data,yerr=data2)
plt.plot(data)
#plt.yscale('log')
plt.show()





