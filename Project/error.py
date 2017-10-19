from numpy import *
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
import statistics 
from array import *
a=[]
b=[]
c=[]
for i in range(0,10):
	a.append([])
	with open('HMC_Results_x_2_'+ str(i)) as inf:
		for line in inf:
			parts = line.split()
			if len(parts) > 0:
				a[i].append(double(parts[0]))

for i in range(0,10):
	print(a[i])

sum=0
sum1=0
std=0

for i in range(0,len(a[0])):

	for j in range(0,len(a)):
		sum += a[j][i]
		sum1 += a[j][i]*a[j][i]
	sum=sum/len(a)
	sum1=sum1/len(a)
	std=(sum1-sum**2)/(len(a)-1)
	b.append(sum)
	c.append(std)
	sum=0
	sum1=0
	std=0

'''
x=np.linspace(0, 2000,num=1999)
plt.errorbar(x,b,yerr=c, linestyle="None",color='blue',fmt='x')
plt.show()
'''

fig = plt.figure()
ax2=fig.add_subplot(211)
line1, = ax2.plot(b,color='red')
ax = fig.add_subplot(212)

line, = ax.plot(c, color='blue', lw=2)

ax.set_yscale('log')
ax2.set_title('Average X^2 for a Harmonic Oscillator')
ax.set_title('Error')
ax2.set_xlabel('Monte Carlo Iterations')
ax.set_xlabel('Monte Carlo Iterations')
ax2.set_ylabel('<x^2>')
ax.set_ylabel('Standard Deviation')

plt.show()




