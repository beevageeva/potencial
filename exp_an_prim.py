from math import pi,e
from scipy import integrate
from exp_common import G
import numpy as np


def potFunc(r, rhoC, r0):
	r[0] = 0.5 * r[1]
	res =  4*G*pi*rhoC*np.exp(-r/r0)* (r0**2 + 2*r0**3/r)
	last = res[res.shape[0] - 1]
	return res - last




def vcFunc(r, rhoC, r0):
	r[0] = 0.5 * r[1]
	#res = 2*np.sqrt(G*pi*rhoC*np.exp(-r/r0)*(-r**2*r0 - 2*r*r0**2 - 2*r0**3)/r)
	res = 2*np.sqrt(G*pi*rhoC*np.exp(-r/r0)*(r**2*r0 + 2*r*r0**2 + 2*r0**3)/r)
	#return res - res[0]
	return res


def massFunc(r, rhoC, r0):
	t1 = e ** (-r/r0)
	res=4*pi*rhoC*np.exp(-r/r0)* (-r**2*r0 - 2*r*r0**2 - 2*r0**3) 
	res = res - res[0]
	last = res[res.shape[0] - 1]
	print "Total mass is %.1e" % last
	return res
