import matplotlib
from scipy import integrate
from math import pi
matplotlib.use('GTKAgg')
import matplotlib.pyplot as plt
from math import sqrt

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
	#G = 1
	#plot function
	
	def f(x):
		#int1 = integrate.quad(lambda y: y**2 * densFunc(y, rhoC, r0), 0, x , limit=NPOINTSQUAD)
		int1 = integrate.quad(lambda y: y**2 * densFunc(y, rhoC, r0), 0, x)
		return (1.0 / x**2) * int1[0]
	
	epsilon = 0.0001
	res = np.zeros(r.shape)
	i = 0
	for x in r:
		#res[i] =  4* pi * G * integrate.quad(f, 0, x, points=[0], limit=NPOINTSQUAD)[0]  ->float division by 0 in 0
		#res[i] =  4* pi * G * integrate.quad(f, epsilon, x, limit=NPOINTSQUAD)[0]
		res[i] =  4* pi * G * integrate.quad(f, epsilon, x)[0]
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
		int1 = integrate.quad(lambda y: y**2 * densFunc(y, rhoC, r0), 0, x)
		t1 =  (4.0 * pi * G * int1[0]) / x

		if(t1<0):
			print("vc**2 for r = %4.2f negative=%4.2f t1=%4.2f, return 0"%(x, pr1, t1))
			res[i] = 0
		else:	
			res[i] = sqrt(t1)
		i+=1
	return res
			 


def dpFunc(r, rhoC, r0):
	res = np.zeros(r.shape)
	i = 0
	for x in r:
		#int1 = integrate.quad(lambda y: (y * densFunc(y, rhoC, r0))/sqrt(y**2-x**2), 0.0001, np.inf)	
		int1 = integrate.quad(lambda y: densFunc(sqrt(y**2 + x**2), rhoC, r0), 0, np.inf)	
		res[i] =  2 *  abs(int1[0]) 
		i+=1
	return res


def massFunc(r, rhoC, r0):
	res = np.zeros(r.shape)
	i = 0
	for x in r:
		int1 = integrate.quad(lambda y: y**2 * densFunc(y, rhoC, r0), 0, x)	
		res[i] =  4.0 * pi *  int1[0]
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
plotForRhoC(rhoCFixed, "dp")


