from math import pi,e
from scipy import integrate
from exp_common import G
import numpy as np


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

			 




def massFunc(r, rhoC, r0):
	t1 = e ** (-r/r0)
	return   4 * pi * rhoC * r0 * (2* r0**2 - 2 * r0 **2 * t1 - 2.0 *r0 * r * t1 - r ** 2 * t1 ) 
	#element by element
#	res = np.zeros(r.shape)
#	i = 0
#	for x in r:
#		t1 = e ** (-x/r0)
#		res[i] =  (-4 * pi * rhoC * r0) * (- 2* r0**2 + 2 * r0 **2 * t1 + 2.0 *r0 * x * t1 + x ** 2 * t1 ) 
#		i+=1
#	return res
