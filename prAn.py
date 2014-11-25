import matplotlib
from scipy import integrate
from math import pi
matplotlib.use('TKAgg')
matplotlib.rcParams['axes.titlesize'] = 'small'
import matplotlib.pyplot as plt
from math import sqrt, e

import numpy as np


def densFunc(r,rhoC, r0):
	return rhoC * np.exp(-r /r0 ) 




def dpFunc(r, rhoC, r0):
	res = np.zeros(r.shape)
	i = 0
	for x in r:
		#int1 = integrate.quad(lambda y: (y * densFunc(y, rhoC, r0))/sqrt(y**2-x**2), 0.0001, np.inf)	
		int1 = integrate.quad(lambda y: densFunc(sqrt(y**2 + x**2), rhoC, r0), 0, np.inf)	
		res[i] =  2 *  abs(int1[0]) 
		i+=1
	return res


def potFunc(r, rhoC, r0):
	#avoid first element of 0
	#eps = 0.0004
	#r[r<eps] = eps
	eps =  0.5 * r[1]
	r[0] = 0.5 * r[1]
	#print("r")
	#print(r)
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
	res=  4 * pi * rhoC * r0 * (2* r0**2 - 2 * r0 **2 * t1 - 2.0 *r0 * r * t1 - r ** 2 * t1 ) 
	last = res[res.shape[0] - 1]
	print "Total mass is %.1e" % last
	return res
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

	plt.ylim(-3e+35,0)	
	#plt.savefig("densNPRM10-%d.png" % int(r0))
	plt.show(block=True)
	#plt.pause(10)

def plotForRhoC(rhoC , r0, fType):
	plt.cla()

	def plotDensFunc():
		plt.plot(r, densFunc(r, rhoC, r0))
		
		plt.xlabel('radius(m)')
		plt.ylabel('density(kg/m3)')
		
		plt.title('Densidad = rhoC * exp(-r/r0) , rhoC = %.1e kg/m3, r0=%.1e m' % (rhoC, r0))


	def plotPotFunc():
		plt.plot(r, potFunc(r, rhoC, r0))
		
		plt.xlabel('radius(m)')
		plt.ylabel('V(J/kg)')
		plt.title('Pot for rho = rhoC * exp(-r/r0), rhoC = %.1e kg/m3, r0=%.1e m' % (rhoC, r0))

	def plotVcFunc():
		plt.plot(r, vcFunc(r, rhoC, r0), label='r0=%.1e' % r0)
		
		plt.xlabel('radius(m)')
		plt.ylabel('vc(m/s)')
		
		plt.title('Vc for rho = rhoC * exp(-r/r0), rhoC = %.1e kg/m3, r0 = %.1e m' % (rhoC, r0))
	
	def plotDpFunc():
		plt.plot(r, dpFunc(r, rhoC, r0), label='r0=%.1e' % r0)
		
		plt.xlabel('radius(m)')
		plt.ylabel('dens. proj.(Kg/m2)')
		
		plt.title('Proj. dens. for rho = rhoC * exp(-r/r0), rhoC = %.1e kg/m3, r0=%.1e m' % (rhoC, r0))

	def plotMassFunc():
		plt.plot(r, massFunc(r, rhoC, r0))
		
		plt.xlabel('radius(m)')
		plt.ylabel('mass(Kg)')
		
		plt.title('Mass for rho = rhoC * exp(-r/r0), rhoC = %.1e kg/m3, r0=%.1e m' % (rhoC, r0))
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



#numPoints = 10000
numPoints = 1000
#Rmax = 10**11
Rmax = 5.0 * 10**11 #for mass
#Rmax = 10**12    #for vc



G = 6.6 * 10 **(-11)
r0Fixed = 1.0*10**10
rhoCFixed = 1.0*10**5



import getopt,sys
try:
				opts, args = getopt.getopt(sys.argv[1:], "", ["r0=",  "rhoC=", "type=", "rmax="])
except getopt.GetoptError as err:
				# print help information and exit:
				print("Error in parsing args")
				print(str(err)) # will print something like "option -a not recognized"
				sys.exit(2)
#print("OPTIONS")
#print(opts)
#print("OPTIONS END")
ptype = "p" 
for o, a in opts:
				#print("o is now >%s<" % o)
				if o == "--type":
								ptype = a
				elif o == "--r0":
								r0Fixed = float(a)
				elif o == "--rhoC":
								rhoCFixed = float(a)
				elif o == "--rmax":
								Rmax = float(a)


				else:
								print("option %s not recognized " % o)

r = np.linspace(0,Rmax,numPoints)
plotForRhoC(rhoCFixed, r0Fixed ,  ptype)


