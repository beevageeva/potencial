import numpy as np
from scipy import integrate
from math import sqrt

G = 6.6 * 10 **(-11)


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
