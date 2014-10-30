import matplotlib
from scipy import integrate
from math import pi
matplotlib.use('TKAgg')
matplotlib.rcParams['axes.titlesize'] = 'small'
import matplotlib.pyplot as plt
from math import sqrt, e

import numpy as np

#numPoints = 10000
numPoints = 1000
#Rmax = 10**11
Rmax = 5.0 * 10**11 #for mass
#Rmax = 10**12    #for vc

#not used, makes no difference!
#NPOINTSQUAD = 100


G = 6.6 * 10 **(-11)
r0Fixed = 1.0*10**10
rhoCFixed = 1.0*10**5



#r = np.arange(0,Rmax,Rmax / numPoints)
r = np.linspace(0,Rmax,numPoints)

def densFunc(r,rhoC, r0):
	return rhoC * np.exp(-r /r0 ) 






def potFunc(r, rhoC, r0):
	eps = 0.0001
	r[r<eps] = eps
	print("r")
	print(r)
	res =  4 * pi * G * rhoC * r0**2 * ((2 * r0)/eps - (e**(-(eps/r0)) * (eps + 2 * r0))/eps + (-2 * r0 + np.exp(-(r/r0))* (r + 2 * r0))/r)
	last = res[res.shape[0] - 1]
	return res - last




def vcFunc(r, rhoC, r0):
   return np.sqrt(4.0 * pi * G  * rhoC *r0*	(2 *  r0**2 / r - np.exp(-(r/r0)) *  (2 * r0**2/r + 2 * r0  + r)))



#	res = np.zeros(r.shape)
#	i = 0
#	epsilon = 0.0001
#	for x in r:
#		if x == 0:
#			x = epsilon
#		t1 =  (-4.0 * pi * G  * rhoC *  r0 ) * (2 * r0 **2 * x **(-2) * e ** (- x/r0) + 2 * r0  * x ** (-1) * e ** (-x/r0) +  e **(-x/r0) - 2 * r0**2/x )	
#
#		if(t1<0):
#			print("vc**2 for r = %4.2f negative=%4.2f t1=%4.2f, return 0"%(x, pr1, t1))
#			res[i] = 0
#		else:	
#			res[i] = sqrt(t1)
#		i+=1
#	return res
			 




def massFunc(r, rhoC, r0):
	t1 = e ** (-r/r0)
	return  4 * pi * rhoC * r0 * (2* r0**2 - 2 * r0 **2 * t1 - 2.0 *r0 * r * t1 - r ** 2 * t1 ) 
	#element by element
#	res = np.zeros(r.shape)
#	i = 0
#	for x in r:
#		t1 = e ** (-x/r0)
#		res[i] =  (-4 * pi * rhoC * r0) * (- 2* r0**2 + 2 * r0 **2 * t1 + 2.0 *r0 * x * t1 + x ** 2 * t1 ) 
#		i+=1
#	return res
 
def plotForAll(plotFunction):
	plotFunction()
	plt.grid(True)
	plt.draw()
	
	#plt.savefig("densNPRM10-%d.png" % int(r0))
	plt.show(block=True)
	#plt.pause(10)

def plotForRhoC(rhoC , r0, fType):
	plt.cla()

	def plotDensFunc():
		plt.plot(r, densFunc(r, rhoC, r0), label='r0=%.1e' % r0)
		
		plt.xlabel('radius(m)')
		plt.ylabel('density(kg/m3)')
		
		plt.title('Densidad = rhoC * exp(-r/r0) , rhoC = %.1e' % rhoC)


	def plotPotFunc():
		plt.plot(r, potFunc(r, rhoC, r0))
		
		plt.xlabel('radius(m)')
		plt.ylabel('V(J/kg)')
		plt.title('Pot for rho = rhoC * exp(-r/r0), rhoC = %.1e, r0=%.1e' % (rhoC, r0))

	def plotVcFunc():
		plt.plot(r, vcFunc(r, rhoC, r0), label='r0=%.1e' % r0)
		
		plt.xlabel('radius(m)')
		plt.ylabel('vc(m/s)')
		
		plt.title('Vc for rho = rhoC * exp(-r/r0), rhoC = %.1e, r0 = %.1e' % (rhoC, r0))
	
	def plotDpFunc():
		plt.plot(r, dpFunc(r, rhoC, r0), label='r0=%.1e' % r0)
		
		plt.xlabel('radius(m)')
		plt.ylabel('dens. proj.(Kg/m2)')
		
		plt.title('Proj. dens. for rho = rhoC * exp(-r/r0) , rhoC = %.1e' % rhoC)

	def plotMassFunc():
		plt.plot(r, massFunc(r, rhoC, r0))
		
		plt.xlabel('radius(m)')
		plt.ylabel('mass(Kg)')
		
		plt.title('Mass for rho = rhoC * exp(-r/r0), rhoC = %.1e, r0=%.1e' % (rhoC, r0))
	if(fType=="p"):
		plotFunc = plotPotFunc
	elif(fType=="v"):
		plotFunc = plotVcFunc
	elif(fType=="m"):
		plotFunc = plotMassFunc
	elif(fType=="d"):
		plotFunc = plotDensFunc
	elif(fType=="dp"):
		plotFunc = plotDpFunc
	else:
		plotFunc = None
		print("undefined function type %s", fType) 
		import sys
		sys.exit(0)
	plotForAll(plotFunc)





plotForRhoC(rhoCFixed, r0Fixed,  "v")


