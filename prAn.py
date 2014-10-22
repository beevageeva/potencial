import matplotlib
from scipy import integrate
from math import pi
matplotlib.use('GTKAgg')
import matplotlib.pyplot as plt
from math import sqrt, e

import numpy as np

#numPoints = 10000
numPoints = 1000
Rmax = 10**11
#Rmax = 5.0 * 10**11 #for mass
#Rmax = 10**12    #for vc

#not used, makes no difference!
#NPOINTSQUAD = 100


G = 6.6 * 10 **(-11)
r0Array = [0.5*10**10, 1.0*10**10, 1.5 * 10**10, 2.0*10**10]
r0Fixed = 1.0*10**10
rhoCArray = [0.5 * 10**5, 1.0*10**5, 1.5*10**5, 2.0*10**5]
rhoCFixed = 1.0*10**5


#legendUp = False
legendUp = True

#r = np.arange(0,Rmax,Rmax / numPoints)
r = np.linspace(0,Rmax,numPoints)

def densFunc(r,rhoC, r0):
	return rhoC * np.exp(-r /r0 ) 



def potFunc(r, rhoC, r0):
	epsilon = 0.0001
	res = np.zeros(r.shape)
	i = 0
	for x in r:
		if x ==0:
			x = epsilon
		int1 = integrate.quad(lambda a: a ** (-2) * e ** (-a/r0) , epsilon, x )[0]
		int2 = integrate.quad(lambda a: a ** (-1) * e ** (-a/r0) , epsilon, x )[0]
		res[i] = (4 * pi * G * rhoC *r0 ** 2)* (2.0 * r0 * (1.0/epsilon - 1.0/x)  -int1 *r0 - 2.0  * int2  + e ** (-x/r0)) 
		i+=1
	return res - res[i-1]
	#not working
	#return 4* pi * G * integrate.quad(f, 0, r)[0]


def vcFunc(r, rhoC, r0):
	res = np.zeros(r.shape)
	i = 0
	epsilon = 0.0001
	for x in r:
		if x == 0:
			x = epsilon
		t1 =  (-4.0 * pi * G  * rhoC *  r0 ) * (2 * r0 **2 * r **(-2) * e ** (- r/r0) + 2 * r0  * r ** (-1) * e ** (-r/r0) +  e **(-r/r0) - 2 * r0**2/r )	

		if(t1<0):
			print("vc**2 for r = %4.2f negative=%4.2f t1=%4.2f, return 0"%(x, pr1, t1))
			res[i] = 0
		else:	
			res[i] = sqrt(t1)
		i+=1
	return res
			 




def massFunc(r, rhoC, r0):
	res = np.zeros(r.shape)
	i = 0
	for x in r:
		t1 = e ** (-B * r)
		res[i] =  (-4 * pi * rhoC * r0) * (- 2* r0**2 + 2 * r0 **2 * t1 + 2.0 *r0 * r * t1 + r ** 2 * t1 ) 
		i+=1
	return res
 
def plotForAll(plotFunction):
	plotFunction()
	plt.grid(True)
	if not legendUp:	
		plt.legend()
	else:
		ax = plt.gca()
		box = ax.get_position()
		ax.set_position([box.x0, box.y0, box.width , box.height*0.85])
		plt.legend(loc='lower center', bbox_to_anchor=(0.5, 1.05),  ncol=2)
	plt.draw()
	
	#plt.savefig("densNPRM10-%d.png" % int(r0))
	plt.show(block=True)
	#plt.pause(10)

def plotForRhoC(rhoC , fType):
	plt.cla()

	def plotDensFunc():
		for r0 in r0Array:
			plt.plot(r, densFunc(r, rhoC, r0), label='r0=%.1e' % r0)
		
		plt.xlabel('radius(m)')
		plt.ylabel('density(kg/m3)')
		
		plt.title('Densidad = rhoC * exp(-r/r0) , rhoC = %.1e' % rhoC)


	def plotPotFunc():
		for r0 in r0Array:
			plt.plot(r, potFunc(r, rhoC, r0), label='r0=%.1e' % r0)
		
		plt.xlabel('radius(m)')
		plt.ylabel('V(J/kg)')
		plt.title('Pot for rho = rhoC * exp(-r/r0) , rhoC = %.1e' % rhoC)

	def plotVcFunc():
		for r0 in r0Array:
			plt.plot(r, vcFunc(r, rhoC, r0), label='r0=%.1e' % r0)
			#print("r0=%.1e,vc(r0)=%.1e" % (r0, vcFunc(np.array([r0]), rhoC, r0)[0]))
		
		plt.xlabel('radius(m)')
		plt.ylabel('vc(m/s)')
		
		plt.title('Vc for rho = rhoC * exp(-r/r0) , rhoC = %.1e' % rhoC)
	
	def plotDpFunc():
		for r0 in r0Array:
			plt.plot(r, dpFunc(r, rhoC, r0), label='r0=%.1e' % r0)
		
		plt.xlabel('radius(m)')
		plt.ylabel('dens. proj.(Kg/m2)')
		
		plt.title('Proj. dens. for rho = rhoC * exp(-r/r0) , rhoC = %.1e' % rhoC)

	def plotMassFunc():
		for r0 in r0Array:
			plt.plot(r, massFunc(r, rhoC, r0), label='r0=%.1e' % r0)
		
		plt.xlabel('radius(m)')
		plt.ylabel('mass(Kg)')
		
		plt.title('Mass for rho = rhoC * exp(-r/r0) , rhoC = %.1e' % rhoC)
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



def plotForRadius(r0, fType):
	
	plt.cla()

	def plotDensFunc():
		for rhoC in rhoCArray:
			plt.plot(r, densFunc(r, rhoC, r0), label='rhoC=%.1e' % rhoC)
		
		plt.xlabel('radius(m)')
		plt.ylabel('density(kg/m3)')
		
		plt.title('Densidad = rhoC * exp(-r/r0) , r0 = %.1e' % r0)


	def plotPotFunc():
		for rhoC in rhoCArray:
			plt.plot(r,  potFunc(r, rhoC, r0), label='rhoC=%.1e' % rhoC)
		
		plt.xlabel('radius(m)')
		plt.ylabel('V(J/kg)')
		plt.title('Pot for rho = rhoC * exp(-r/r0) , r0 = %.1e' % r0)

	def plotVcFunc():
		for rhoC in rhoCArray:
			plt.plot(r, vcFunc(r, rhoC, r0), label='rhoC=%.1e' % rhoC)
		#print("r0=%.1e,vc(r0)=%.1e" % (r0, vcFunc(np.array([r0]), rhoC, r0)[0]))
		plt.xlabel('radius(m)')
		plt.ylabel('vc(m/s)')
		
		plt.title('Vc for rho = rhoC * exp(-r/r0) , r0 = %.1e' % r0)

	def plotDpFunc():
		for rhoC in rhoCArray:
			plt.plot(r, dpFunc(r, rhoC, r0), label='rhoC=%.1e' % rhoC)
		
		plt.xlabel('radius(m)')
		plt.ylabel('Proj dens(kg/m2)')
		
		plt.title('Proj. dens. for rho = rhoC * exp(-r/r0) , r0 = %.1e' % r0)
	
	def plotMassFunc():
		for rhoC in rhoCArray:
			plt.plot(r, massFunc(r, rhoC, r0), label='rhoC=%.1e' % rhoC)
		
		plt.xlabel('radius(m)')
		plt.ylabel('Mass(Kg)')
		
		plt.title('Mass for rho = rhoC * exp(-r/r0) , r0 = %.1e' % r0)
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



import datetime as dt

#plotForRadius(r0Fixed, "dp")
plotForRhoC(rhoCFixed, "p")


