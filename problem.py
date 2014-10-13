import matplotlib
from scipy import integrate
from math import pi
matplotlib.use('GTKAgg')
import matplotlib.pyplot as plt
from math import sqrt

import numpy as np

numPoints = 1000
domMult = 10000


def densFunc(r,rhoC, r0):
	return rhoC * np.exp(-r /r0 ) 



def potFunc(r, rhoC, r0):
	G = 6.6 * 10 **(-11)
	#G = 1
	#plot function
	
	def f(x):
		int1 = integrate.quad(lambda y: y**2 * densFunc(y, rhoC, r0), 0, x )
		return (1.0 / x**2) * int1[0]
	

	res = np.zeros(r.shape)
	i = 0
	for x in r:
		res[i] =  4* pi * G * integrate.quad(f, 0.0001, x)[0]
		i+=1
	return res
	#not working
	#return 4* pi * G * integrate.quad(f, 0, r)[0]


def vcFunc(r, rhoC, r0):
	G = 6.6 * 10 **(-11)
	res = np.zeros(r.shape)
	i = 0
	for x in r:
		if x == 0:
			x = 0.0001
		int1 = integrate.quad(lambda y: y**2 * densFunc(y, rhoC, r0), 0, x )
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
		int1 = integrate.quad(lambda y: (y * densFunc(y, rhoC, r0))/sqrt(y**2-x**2), x, np.inf )	
		res[i] =  2 *  int1[0] 
		i+=1
	return res


def massFunc(r, rhoC, r0):
	res = np.zeros(r.shape)
	i = 0
	for x in r:
		int1 = integrate.quad(lambda y: y**2 * densFunc(y, rhoC, r0), 0, x )	
		res[i] =  4.0 * pi *  int1[0]
		i+=1
	return res
 


def plotForRhoC(rhoC):
	#r0Array = [0.5, 1.0, 2.0]
	r0Array = [1.0]
	maxX = max(r0Array) * domMult
	r = np.arange(0,maxX,maxX / numPoints)
	plt.cla()

	def plotDensFunc():
		for r0 in r0Array:
			plt.plot(r, densFunc(r, rhoC, r0), label='r0=%4.1f' % r0)
		
		plt.xlabel('radius')
		plt.ylabel('density')
		
		plt.title('Densidad = rhoC * exp(-r/r0) , rhoC = %5.1f' % rhoC)


	def plotPotFunc():
		for r0 in r0Array:
			plt.plot(r, potFunc(r, rhoC, r0), label='r0=%5.1f' % r0)
		
		plt.xlabel('radius')
		plt.ylabel('pot')
		
		plt.title('Pot for rho = rhoC * exp(-r/r0) , rhoC = %5.1f' % rhoC)

	def plotVcFunc():
		for r0 in r0Array:
			plt.plot(r, vcFunc(r, rhoC, r0), label='r0=%5.1f' % r0)
		
		plt.xlabel('radius')
		plt.ylabel('vc')
		
		plt.title('Vc for rho = rhoC * exp(-r/r0) , rhoC = %5.1f' % rhoC)
	
	def plotDpFunc():
		for r0 in r0Array:
			plt.plot(r, dpFunc(r, rhoC, r0), label='r0=%5.1f' % r0)
		
		plt.xlabel('radius')
		plt.ylabel('dens. proj.')
		
		plt.title('Proj. dens. for rho = rhoC * exp(-r/r0) , rhoC = %5.1f' % rhoC)

	def plotMassFunc():
		for r0 in r0Array:
			plt.plot(r, massFunc(r, rhoC, r0), label='r0=%5.1f' % r0)
		
		plt.xlabel('radius')
		plt.ylabel('mass')
		
		plt.title('Mass for rho = rhoC * exp(-r/r0) , rhoC = %5.1f' % rhoC)
	#plotDensFunc()
	plotPotFunc()
	#plotVcFunc()
	#plotDpFunc()
	#plotMassFunc()
	plt.grid(True)
	plt.legend()
	plt.draw()
	
	#plt.savefig("densNPRM10-%d.png" % int(r0))
	plt.show(block=True)
	#plt.pause(10)



def plotForRadius(r0):
	
	rhoCArray = [0.5, 1.0, 2.0]
	
	r = np.arange(0,r0*domMult,r0*domMult / numPoints)
	plt.cla()

	def plotDensFunc():
		for rhoC in rhoCArray:
			plt.plot(r, densFunc(r, rhoC, r0), label='rhoC=%2.1f' % rhoC)
		
		plt.xlabel('radius')
		plt.ylabel('density')
		
		plt.title('Densidad = rhoC * exp(-r/r0) , r0 = %5.1f' % r0)


	def plotPotFunc(plotFirst, plotSecond):
		for rhoC in rhoCArray:
			plt.plot(r, potFunc(r, rhoC, r0), label='rhoC=%2.1f' % rhoC)
		
		plt.xlabel('radius')
		plt.ylabel('pot')
		
		plt.title('Pot for rho = rhoC * exp(-r/r0) , r0 = %5.1f' % r0)

	def plotVcFunc():
		for rhoC in rhoCArray:
			plt.plot(r, vcFunc(r, rhoC, r0), label='rhoC=%2.1f' % rhoC)
		plt.xlabel('radius')
		plt.ylabel('vc')
		
		plt.title('Vc for rho = rhoC * exp(-r/r0) , r0 = %5.1f' % r0)

	def plotDpFunc():
		for rhoC in rhoCArray:
			plt.plot(r, dpFunc(r, rhoC, r0), label='rhoC=%2.1f' % rhoC)
		
		plt.xlabel('radius')
		plt.ylabel('Proj dens')
		
		plt.title('Proj. dens. for rho = rhoC * exp(-r/r0) , r0 = %5.1f' % r0)
	
	def plotMassFunc():
		for rhoC in rhoCArray:
			plt.plot(r, massFunc(r, rhoC, r0), label='rhoC=%2.1f' % rhoC)
		
		plt.xlabel('radius')
		plt.ylabel('Mass')
		
		plt.title('Mass for rho = rhoC * exp(-r/r0) , r0 = %5.1f' % r0)
	#plotDensFunc()
	#plotPotFunc(True, True)
	#plotPotFunc(False, True)
	#plotPotFunc(False, False)
	#plotVcFunc()
	#plotDpFunc()
	plotMassFunc()
	plt.grid(True)
	plt.legend()
	plt.draw()
	
	#plt.savefig("densNPRM10-%d.png" % int(r0))
	plt.show(block=True)
	#plt.pause(10)


def plotDensFuncDifR():
	plt.title('rho = rhoC * exp(-r/r0) , rhoC = 1')
	for r0 in (1.0, 100.0, 10000.0):
		plotForRadius(r0)



#plotForRadius(1.0)
#plotForRadius(100.0)
#plotForRadius(10000.0)
plotForRhoC(1.0)


#for r0 in (1.0, 100.0, 10000.0):
#	plotForRadius(r0)
