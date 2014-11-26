import matplotlib
from scipy import integrate
from math import pi
#if the following is causing Segmentation fault see bug:
#http://sourceforge.net/p/matplotlib/mailman/message/12934184/
#matplotlib.use('GTKAgg')
#use
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
from math import sqrt

import numpy as np


from exp_common import densFunc, dpFunc
#to use numerical functions
#from exp_num import vcFunc, potFunc, massFunc
#to use analytical functions
from exp_an import vcFunc, potFunc, massFunc



 
def plotForAll(plotFunction):
	plt.cla()
	plotFunction()
	plt.grid(True)
	#to set ylim
	#plt.ylim(-3.5e+35,0)
	#plt.ylim(0,200000)
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

def plotForRhoC(rhoC, r0Array,fType):

	def plotDensFunc():
		for r0 in r0Array:
			plt.plot(r, densFunc(r, rhoC, r0), label='r0=%.1e m' % r0)
		
		plt.xlabel('radius(m)')
		plt.ylabel('density(kg/m3)')
		
		plt.title('Densidad = rhoC * exp(-r/r0) , rhoC = %.1e kg/m3' % rhoC)


	def plotPotFunc():
		for r0 in r0Array:
			plt.plot(r, potFunc(r, rhoC, r0), label='r0=%.1e m' % r0)
		
		plt.xlabel('radius(m)')
		plt.ylabel('V(J/kg)')
		plt.title('Pot for rho = rhoC * exp(-r/r0) , rhoC = %.1e kg/m3' % rhoC)

	def plotVcFunc():
		for r0 in r0Array:
			plt.plot(r, vcFunc(r, rhoC, r0), label='r0=%.1e m' % r0)
			#print("r0=%.1e,vc(r0)=%.1e" % (r0, vcFunc(np.array([r0]), rhoC, r0)[0]))
		
		plt.xlabel('radius(m)')
		plt.ylabel('vc(m/s)')
		
		plt.title('Vc for rho = rhoC * exp(-r/r0) , rhoC = %.1e kg/m3' % rhoC)
	
	def plotDpFunc():
		for r0 in r0Array:
			plt.plot(r, dpFunc(r, rhoC, r0), label='r0=%.1e m' % r0)
		
		plt.xlabel('radius(m)')
		plt.ylabel('dens. proj.(Kg/m2)')
		
		plt.title('Proj. dens. for rho = rhoC * exp(-r/r0) , rhoC = %.1e kg/m3' % rhoC)

	def plotMassFunc():
		for r0 in r0Array:
			plt.plot(r, massFunc(r, rhoC, r0), label='r0=%.1e m' % r0)
		
		plt.xlabel('radius(m)')
		plt.ylabel('mass(Kg)')
		
		plt.title('Mass for rho = rhoC * exp(-r/r0) , rhoC = %.1e kg/m3' % rhoC)
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



def plotForRadius(r0, rhoCArray, fType):
	

	def plotDensFunc():
		for rhoC in rhoCArray:
			plt.plot(r, densFunc(r, rhoC, r0), label='rhoC=%.1e kg/m3' % rhoC)
		
		plt.xlabel('radius(m)')
		plt.ylabel('density(kg/m3)')
		
		plt.title('Densidad = rhoC * exp(-r/r0) , r0 = %.1e m' % r0)


	def plotPotFunc():
		for rhoC in rhoCArray:
			plt.plot(r,  potFunc(r, rhoC, r0), label='rhoC=%.1e kg/m3' % rhoC)
		
		plt.xlabel('radius(m)')
		plt.ylabel('V(J/kg)')
		plt.title('Pot for rho = rhoC * exp(-r/r0) , r0 = %.1e m' % r0)

	def plotVcFunc():
		for rhoC in rhoCArray:
			plt.plot(r, vcFunc(r, rhoC, r0), label='rhoC=%.1e kg/m3' % rhoC)
		#print("r0=%.1e,vc(r0)=%.1e" % (r0, vcFunc(np.array([r0]), rhoC, r0)[0]))
		plt.xlabel('radius(m)')
		plt.ylabel('vc(m/s)')
		
		plt.title('Vc for rho = rhoC * exp(-r/r0) , r0 = %.1e m' % r0)

	def plotDpFunc():
		for rhoC in rhoCArray:
			plt.plot(r, dpFunc(r, rhoC, r0), label='rhoC=%.1e kg/m3' % rhoC)
		
		plt.xlabel('radius(m)')
		plt.ylabel('Proj dens(kg/m2)')
		
		plt.title('Proj. dens. for rho = rhoC * exp(-r/r0) , r0 = %.1e m' % r0)
	
	def plotMassFunc():
		for rhoC in rhoCArray:
			plt.plot(r, massFunc(r, rhoC, r0), label='rhoC=%.1e kg/m3' % rhoC)
		
		plt.xlabel('radius(m)')
		plt.ylabel('Mass(Kg)')
		
		plt.title('Mass for rho = rhoC * exp(-r/r0) , r0 = %.1e m' % r0)
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
r0Array = [0.5*10**10, 1.0*10**10, 1.5 * 10**10, 2.0*10**10]
r0Fixed = 1.0*10**10
rhoCArray = [0.5 * 10**5, 1.0*10**5, 1.5*10**5, 2.0*10**5]
rhoCFixed = 1.0*10**5

#if true it will plot the legend in a box in the upper part
#legendUp = False
legendUp = True


fixedType = "r0" #or rhoC


import getopt,sys
try:
				opts, args = getopt.getopt(sys.argv[1:], "", ["r0=",  "rhoC=", "type=", "rmax=", "fixed="])
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
								arg1 = a.split(",")
								if(len(arg1)==1):
									r0Fixed = float(a)
								else:
									r0Array = []
									for a1 in arg1:
										r0Array.append(float(a1))
				elif o == "--rhoC":
								arg1 = a.split(",")
								if(len(arg1)==1):
									rhoCFixed = float(a)
								else:
									rhoCArray = []
									for a1 in arg1:
										rhoCArray.append(float(a1))
				elif o == "--rmax":
								Rmax = float(a)
				elif o == "--fixed":
								fixedType = a


				else:
								print("option %s not recognized " % o)



#r = np.arange(0,Rmax,Rmax / numPoints)
r = np.linspace(0,Rmax,numPoints)

if fixedType == "r0":
	plotForRadius(r0Fixed, rhoCArray,  ptype)
elif fixedType == "rhoC":
	plotForRhoC(rhoCFixed, r0Array, ptype)
else:
	print("Unknown fixed type %s" % fixedType)



