from numpy import *
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
import statistics 
from array import *

avgx=[];avgx2=[];action=[];KE=[];delta_h=[];i=[];
Mavgx=[];Mavgx2=[];Maction=[];MKE=[];Mdelta_h=[];
avgxerr=[];avgx2err=[];actionerr=[];KEerr=[];delta_herr=[];
sum1=0;sum2=0;sum3=0;sum4=0;sum5=0
sum12=0;sum22=0;sum32=0;sum42=0;sum52=0


file = open("HMC_Stats.dat",'r')

for line in file:
	a,b,c,d,e,f,gg,hh = line.split(' ', 7)
	i.append(float(a))
	avgx.append(float(b))
	avgxerr.append(float(c))
	avgx2.append(float(d))
	avgx2err.append(float(e))
	action.append(float(f))
	KE.append(float(gg))
	delta_h.append(float(hh))

for j in range(1,len(i)):
	sum1 += avgx[j]
	sum2 += avgx2[j]
	sum3 += action[j]
	sum4 += KE[j]
	sum5 += delta_h[j]
	
	sum12 += avgx[j] * avgx[j]
	sum22 += avgx2[j] * avgx2[j]
	sum32 += action[j] * action[j]
	sum42 += KE[j] * KE[j]
	sum52 += delta_h[j] * delta_h[j]
	
	Mavgx.append(sum1/j)
	Mavgx2.append(sum2/j)
	Maction.append(sum3/j)
	MKE.append(sum4/j)
	Mdelta_h.append(sum5/j)

	avgxerr.append((sum12/j)-(sum1*sum1/(j*j))/j)
	avgx2err.append((sum22/j)-(sum2*sum2/(j*j))/j)
	actionerr.append((sum32/j)-(sum3*sum3/(j*j))/j)
	KEerr.append((sum42/j)-(sum4*sum4/(j*j))/j)
	delta_herr.append((sum52/j)-(sum5*sum5/(j*j))/j)






g=plt.figure()
plt.ylabel('<X>')
plt.xlabel('Monte Carlo Iterations')
plt.title('Average X')
plt.plot(Mavgx,label='Simulated AvgX')
plt.axhline(y=0,lw=1,color='red',label='Theoritical AvgX')
plt.legend(loc='upper right')
g.savefig("Average_X.pdf")


h=plt.figure()
plt.ylabel('<X^2>')
plt.xlabel('Monte Carlo Iterations')
plt.title('Average X^2')
hx2=h.add_subplot(211)
hx2.plot(Mavgx2,label='Simulated AvgX^2')
hx2.axhline(y=0.4472135955,lw=1,color='red',label='Theoritical AvgX^2')
hx2.legend(loc='upper right')
#hx=h.add_subplot(212)
#hx.plot(avgx2err,label='Simulated AvgX^2 Error')
h.savefig("Average_X2.pdf")
print(avgx2err)

k=plt.figure()
plt.xlabel('Monte Carlo Iterations')
plt.ylabel('Action Per Lattice Site(S)')
plt.title('Average Action')
plt.plot(Maction,label='Simulated Action')
plt.axhline(y=0.5,lw=1,color='red',label='Theoritical Action')
plt.legend(loc='upper right')
k.savefig("Average_Action.pdf")


l=plt.figure()
plt.xlabel('Monte Carlo Iterations')
plt.ylabel('Kinetic Energy Per Lattice Site')
plt.title('Average Kinetic Energy')
plt.plot(MKE,label='Simulated Avgerage Kinetic Energy')
plt.axhline(y=0.5,lw=1,color='red',label='Theoritical Avgerage Kinetic Energy')
plt.legend(loc='lower right')
l.savefig("Average_KE.pdf")

m=plt.figure()
plt.xlabel('Monte Carlo Iterations')
plt.ylabel('Delta_H')
plt.title('Average Delta_H')
plt.plot(Mdelta_h,label='Simulated Delta_H')
plt.axhline(y=0,lw=1,color='red',label='Theoritical Delta_H')
plt.legend(loc='upper right')
m.savefig("Average_Delta_H.pdf")

n=plt.figure()
x = np.linspace(-2,2,100) # 100 linearly spaced numbers
y = (1/(3.141**0.5))*np.exp(-x**2)
#out1 = plt.plot(x,y)
out = plt.hist(avgx,bins=100,normed=1)
n.savefig("Wavefunction.pdf")














