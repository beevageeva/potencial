import matplotlib
matplotlib.use('TKAgg')
from scipy import integrate
from math import pi
import matplotlib.pyplot as plt
from math import sqrt
import numpy as np

def densFunc(r,rhoC, r0):
	return rhoC * np.exp(-r /r0 ) 

epsilon = 0.0004
#epsilon = 1



def func1Arg(r, rhoC, r0):
	res = np.zeros(r.shape)
	i = 0
	for x in r:
		if(x==0):
			res[i] = 0
		else:
			res[i] = x**(-1) * densFunc(x, rhoC, r0)
		i+=1
	return res

def func2Arg(r, rhoC, r0):
	res = np.zeros(r.shape)
	i = 0
	for x in r:
		if(x==0):
			res[i] = 0
		else:
			res[i] = x**(-2) * densFunc(x, rhoC, r0)
		i+=1
	return res

def func1(r, rhoC, r0):
	def f(x):
		if x==0:
			return 0
		int1 = integrate.quad(lambda y: y**(-1) * densFunc(y, rhoC, r0), epsilon, x )
		return int1[0]
	res = np.zeros(r.shape)
	i = 0
	for x in r:
		res[i] = f(x)
		i+=1
	return res
	#not working
	#return f(r)

def func2(r, rhoC, r0):
	def f(x):
		if x==0:
			return 0
		int1 = integrate.quad(lambda y: y**(-2) * densFunc(y, rhoC, r0), epsilon, x )
		return int1[0]
	res = np.zeros(r.shape)
	i = 0
	for x in r:
		res[i] = f(x)
		i+=1
	return res
	#not working
	#return f(r)
	
def anFunc1I2(r, rhoC, r0):
	t1 = np.exp(-r/r0)
	return -rhoC * ( 2 * r0 **3 * t1 + 2 * r0 **2 * r * t1 + r0 * r **2  * t1  - 2* r0 **3)


def plotFunc():
	domMult = 2
	numPoints = 1000 
	r0 = 0.5
	rhoC = 1.0
	r = np.arange(0,r0*domMult,r0*domMult / numPoints)
	#function 1 
	#plt.title("r**2 * rhoC * exp(-r/r0) , rhoC = %2.1f, r0 = %2.1f" % (rhoC, r0))
	#plt.plot(r, func1(r, rhoC, r0))
	#function int
	#plt.title("Integral(0, r, x**(-2) * rhoC * exp(-x/r0)) , rhoC = %2.1f, r0 = %2.1f" % (rhoC, r0))
	#plt.title("Integral(0, r, x**(-1) * rhoC * exp(-x/r0)) , rhoC = %2.1f, r0 = %2.1f" % (rhoC, r0))
	#plt.title(" x**(-2) * rhoC * exp(-x/r0) , rhoC = %2.1f, r0 = %2.1f" % (rhoC, r0))
	plt.title(" x**(-1) * rhoC * exp(-x/r0) , rhoC = %2.1f, r0 = %2.1f" % (rhoC, r0))
	#plt.plot(r, func2Arg(r, rhoC, r0))
	plt.plot(r, func1Arg(r, rhoC, r0))
	#plt.plot(r, func1(r, rhoC, r0))




plt.clf()
plotFunc()

plt.show(block=True)

