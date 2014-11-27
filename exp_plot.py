import matplotlib
from math import pi
matplotlib.use('TKAgg')
matplotlib.rcParams['axes.titlesize'] = 'small'
import matplotlib.pyplot as plt
from math import sqrt, e

import numpy as np





 
def plotForAll(plotFunction):
	plotFunction()
	plt.grid(True)
	plt.draw()
	#YLIM use in potential case
	#plt.ylim(-3e+35,0)	
	#plt.savefig("densNPRM10-%d.png" % int(r0))
	plt.show(block=True)
	#plt.pause(10)

def plotFor(rhoC , r0, fType):
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



r0Fixed = 1.0*10**10
rhoCFixed = 1.0*10**5


numericalSol = False



import getopt,sys
try:
				opts, args = getopt.getopt(sys.argv[1:], "", ["r0=",  "rhoC=", "type=", "rmax=", "numerical"])
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
				elif o == "--numerical":
								numericalSol = True
				else:
								print("option %s not recognized " % o)

from exp_common import densFunc, dpFunc
if numericalSol:
	#to use numerical functions
	print("using numerical solutions for mass, velocity and potential")	
	from exp_num import vcFunc, potFunc, massFunc
else:
	#to use analytical functions
	print("using analytical solutions for mass, velocity and potential")	
	from exp_an import vcFunc, potFunc, massFunc
	#from exp_an_prim import vcFunc, potFunc, massFunc


r = np.linspace(0,Rmax,numPoints)
plotFor(rhoCFixed, r0Fixed ,  ptype)


