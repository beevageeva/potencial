import matplotlib
from scipy import integrate
from math import pi
matplotlib.use('TKAgg')
matplotlib.rcParams['axes.titlesize'] = 'small'
import matplotlib.pyplot as plt
from math import sqrt, e

import numpy as np


def densFunc(r,b, M):
	a = np.sqrt(b**2 + r**2)
	return M * ((3*(b+a)*a**2 - r**2*(b+3*a) )/(4*pi*(b+a)**3 * a**3 ))	
	


def dpFunc(r, b, M):
	res = np.zeros(r.shape)
	i = 0
	for x in r:
		int1 = integrate.quad(lambda y: densFunc(sqrt(y**2 + x**2), b, M), 0, np.inf)	
		res[i] =  2 *  abs(int1[0]) 
		i+=1
	return res


def potFunc(r, b, M):
	return -G * M / (b + np.sqrt(b**2 + r**2))

def vcFunc(r, b, M):
	a = np.sqrt(b**2 + r**2)	
	return np.sqrt((G * M * r ** 2)/((b + a)**2 * a))


def massFunc(r, b, M):
	return M*r**3/((b + np.sqrt(b**2 + r**2))**2*np.sqrt(b**2 + r**2))
 
def plotForAll(plotFunction):
	plotFunction()
	plt.grid(True)
	plt.draw()
	#YLIM
	#plt.ylim(0,100000)
	#plt.ylim(0,9e-21)
	
	#plt.savefig("densNPRM10-%d.png" % int(r0))
	plt.show(block=True)
	#plt.pause(10)

def plotFor(b , M, fType):
	plt.cla()

	def plotDensFunc():
		plt.plot(r, densFunc(r, b, M))
		
		plt.xlabel('radius(m)')
		plt.ylabel('density(kg/m3)')
		
		plt.title('Densidad  b=%.1e m, M = %.1e kg' % (b, M))


	def plotPotFunc():
		plt.plot(r, potFunc(r, b, M))
		
		plt.xlabel('radius(m)')
		plt.ylabel('V(J/kg)')
		plt.title('Pot b = %.1e m, M=%.1e Kg' % (b, M))

	def plotVcFunc():
		vcRes = vcFunc(r, b, M)
		plt.plot(r, vcRes)
		#mark points for sun velocity
		sunDistance = 2.5e20
		indexR = np.argmin(np.abs(sunDistance - r))	
		print("vc at sunDistance %.2e = %.2e" % (r[indexR], vcRes[indexR]))
		#YLIM
		print("Using ylim")
		plt.ylim(0,300000)	
		#plt.xticks(list(plt.xticks()[0]) +[r[indexR]])
		#remove first tick
		newxticks = list(plt.xticks()[0])
		newxticks.pop(0)
		plt.xticks(newxticks +[r[indexR]])
		plt.yticks(list(plt.yticks()[0]) + [vcRes[indexR]])

	
		plt.xlabel('radius(m)')
		plt.ylabel('vc(m/s)')
		
		plt.title('Vc b = %.1e m M = %.1e kg' % (b, M))
	
	def plotDpFunc():
		plt.plot(r, dpFunc(r, b, M))
		
		plt.xlabel('radius(m)')
		plt.ylabel('dens. proj.(Kg/m2)')
		
		plt.title('Proj. dens. b = %.1e m, M=%.1e Kg' % (b, M))

	def plotMassFunc():
		plt.plot(r, massFunc(r, b, M))
		
		plt.xlabel('radius(m)')
		plt.ylabel('mass(Kg)')
		
		plt.title('Mass b = %.1e m, M=%.1e kg' % (b, M))
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
#bFixed = 1.0*10**10
MFixed = 1.0*10**5



import getopt,sys
try:
				#opts, args = getopt.getopt(sys.argv[1:], "", ["b=",  "rhoC=", "type=", "rmax="])
				opts, args = getopt.getopt(sys.argv[1:], "", ["mass=",  "rhoC=", "type=", "rmax="])
except getopt.GetoptError as err:
				# print help information and exit:
				print("Error in parsing args")
				print(str(err)) # will print something like "option -a not recognized"
				sys.exit(2)
#print("OPTIONS")
#print(opts)
#print("OPTIONS END")
ptype = "p"
rhoC=1.0*10**5 
for o, a in opts:
				#print("o is now >%s<" % o)
				if o == "--type":
								ptype = a
				elif o == "--mass":
								MFixed = float(a)
#				elif o == "--b":
#								bFixed = float(a)
				elif o == "--rhoC":
								rhoC = float(a)
				elif o == "--rmax":
								Rmax = float(a)


				else:
								print("option %s not recognized " % o)

#MFixed = 16.0 * pi * bFixed ** 3 * rhoC / 3.0
bFixed = ((MFixed * 3.0) / (16.0 * pi * rhoC))**(1.0/3)
r = np.linspace(0,Rmax,numPoints)
print("Total mass = %.1e , rhoC = %.1e, b = %.1e" % (MFixed, rhoC, bFixed))
plotFor(bFixed, MFixed ,  ptype)


