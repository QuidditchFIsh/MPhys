from numpy import *
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
import statistics 
from array import *
import math
from matplotlib import gridspec

#initalise arrays and variables
avgx=[];avgx2=[];action=[];KE=[];delta_h=[];i=[];avgx4=[];jjj=[]
Mavgx=[];Mavgx2=[];Maction=[];MKE=[];Mdelta_h=[];Mavgx4=[]
avgxerr=[];avgx2err=[];actionerr=[];KEerr=[];delta_herr=[];avgx4err=[]
sum1=0;sum2=0;sum3=0;sum4=0;sum5=0;sum6=0
sum12=0;sum22=0;sum32=0;sum42=0;sum52=0;sum62=0

#import all data from the file
file  = open("HMC_Stats.dat",'r')
#file1 = open("HMC_X.dat","r")

for line in file:
	a,b,c,d,e,f,gg = line.split(' ', 6)
	i.append(float(a))
	avgx.append(float(b))
	#avgxerr.append(float(c))
	avgx2.append(float(d))
	avgx2err.append(float(e))
	action.append(float(f))
	KE.append(float(gg))
	#delta_h.append(float(hh))
	#avgx4.append(float(ii))

data=np.genfromtxt("HMC_X.dat", unpack=True)
#print(data)


po=[]
#stats calculations
for j in range(0,len(i)):
	sum1 += avgx[j]
	sum2 += avgx2[j]
	sum3 += action[j]
	sum4 += KE[j]
	#sum5 += delta_h[j]
	
	sum12 += avgx[j] * avgx[j]
	sum22 += avgx2[j] * avgx2[j]
	sum32 += action[j] * action[j]
	sum42 += KE[j] * KE[j]
	#sum52 += delta_h[j] * delta_h[j]

	k=j+1
	Mavgx.append(sum1/k)
	Mavgx2.append(sum2/k)
	Maction.append(sum3/k)
	MKE.append(sum4/k)
	#Mdelta_h.append(sum5/j)


	avgxerr.append((sum12/k)-(sum1*sum1/(k*k))/k)
	#avgx2err.append((sum22/k)-(sum2*sum2/(k*k))/k)
	actionerr.append((sum32/k)-(sum3*sum3/(k*k))/k)
	KEerr.append((sum42/k)-(sum4*sum4/(k*k))/k)
	#delta_herr.append((sum52/j)-(sum5*sum5/(j*j))/j)

print(Mavgx[-1])
print(Mavgx2[-1])



#Plotting
g=plt.figure(figsize=(8,6))
gs1 = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
gx2=plt.subplot(gs1[0])
plt.title('Average X')
plt.ylabel('<X>')
gx2.axhline(y=0.0,lw=1,color='red',label='Theoritical AvgX')
gx2.plot(Mavgx,label='Simulated AvgX')
gx2.set_ylim([-0.1,0.1])
plt.legend(loc='lower right')

gx=plt.subplot(gs1[1])
plt.ylabel('<X> error')
plt.xlabel('Iterations')
gx.plot(avgxerr,label='Simulated AvgX Error',color='green')
g.savefig("Average_X_Anharmonic.pdf")



h=plt.figure(figsize=(8,6))
gs2 = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
hx2=plt.subplot(gs2[0])
plt.title('Average X^2')
plt.ylabel('<X^2>')
#hx2.axhline(y=0.4472135955,lw=1,color='red',label='Theoritical AvgX^2')
hx2.plot(Mavgx2,label='Simulated AvgX^2')
#hx2.set_ylim([0,0.5])
plt.legend(loc='center right')

hx=plt.subplot(gs2[1])
plt.ylabel('<X^2> error')
plt.xlabel('Iterations')
hx.plot(avgx2err,label='Simulated AvgX^2 Error',color='green')
h.savefig("Average_X2_Anharmonic.pdf")


o=plt.figure()
plt.ylabel('<X^4>')
plt.xlabel('Monte Carlo Iterations')
plt.title('Average X^4')
#hx2=h.add_subplot(211)
plt.plot(po,label='Simulated AvgX^4')
#plt.axhline(y=0.4472135955,lw=1,color='red',label='Theoritical AvgX^2')
plt.legend(loc='center right')
#hx=h.add_subplot(212)
#hx.plot(avgx2err,label='Simulated AvgX^2 Error')
o.savefig("Average_X4_Anharmonic.pdf")
#o.savefig("Average_X2_harmonic.pdf")

k=plt.figure()
axes=plt.gca()
plt.xlabel('Monte Carlo Iterations')
plt.ylabel('Action Per Lattice Site(S)')
plt.title('Average Action')
plt.plot(Maction,label='Simulated Action')
plt.axhline(y=0.5,lw=1,color='red',label='Theoritical Action')
plt.legend(loc='lower right')
axes.set_ylim([0,0.6])
k.savefig("Average_Action_Anharmonic.pdf")
#k.savefig("graphs/Average_Action_harmonic.pdf")


l=plt.figure()
plt.xlabel('Monte Carlo Iterations')
plt.ylabel('Kinetic Energy Per Lattice Site')
plt.title('Average Kinetic Energy')
plt.plot(MKE,label='Simulated Avgerage Kinetic Energy')
plt.axhline(y=0.5,lw=1,color='red',label='Theoritical Avgerage Kinetic Energy')
plt.legend(loc='upper right')
l.savefig("Average_KE_Anharmonic.pdf")
#l.savefig("graphs/Average_KE_harmonic.pdf")
'''
m=plt.figure()
plt.xlabel('Monte Carlo Iterations')
plt.ylabel('Delta_H')
plt.title('Average Delta_H')
plt.plot(Mdelta_h,label='Simulated Delta_H')
plt.axhline(y=0,lw=1,color='red',label='Theoritical Delta_H')
plt.legend(loc='lower right')
m.savefig("Average_Delta_H_Anharmonic.pdf")
#m.savefig("graphs/Average_Delta_H_harmonic.pdf")
'''

n=plt.figure()
#x = np.linspace(-2,2,100) # 100 linearly spaced numbers
#y = (1/(3.141**0.5))*np.exp(-x**2)
#out1 = plt.plot(x,y)
#out = plt.hist(data,bins=150,normed =1)
n.savefig("Wavefunction_Anharmonic.pdf")
#n.savefig("graphs/Wavefunction_harmonic.pdf")













