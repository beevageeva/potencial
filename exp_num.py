import numpy as np
from math import pi, sqrt
from scipy import integrate


from exp_common import densFunc, G

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
			 




def massFunc(r, rhoC, r0):
	res = np.zeros(r.shape)
	i = 0
	for x in r:
		int1 = integrate.quad(lambda y: y**2 * densFunc(y, rhoC, r0), 0, x)	
		res[i] =  4.0 * pi *  int1[0]
		i+=1
	return res
