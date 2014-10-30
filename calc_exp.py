from scipy import integrate
from math import e, sqrt, pi, exp
from numpy import inf
import sys
#for the potential there will be two integration constants K[0] and K[1]
#for the rest(velocity , masss, integrated mass) there is only one K[0]
#numvars is the hash containing variable parameters(A, B, ...)

G = 6.6 * 10 **(-11)

def testNumvars(numvars):
	if(not 'A' in numvars.keys() or not 'B' in numvars.keys() ):
		print("parameters A and B not present in numvars:")
		print(numvars)
		sys.exit(0)
	


# r here is only one element of radius array used to plot
def calculateP(numvars,K,r):
	epsilon = 0.0004
	if(r==0):
		r = epsilon
	testNumvars(numvars)
	B = numvars["B"]
	A = numvars["A"]
	if r<epsilon:
		r = epsilon
	result = (4 * pi * G * A / B ** 2)* (2.0/(B * epsilon) - (e**(-( B * epsilon)) * (epsilon + 2.0 /B ))/epsilon + (-2.0/B + exp(-(r*B))* (r + 2.0/B))/r) + K[0] / r  + K[1]
	return result

def calculateV(numvars, K , r):
	if(r==0):
		return 0
	testNumvars(numvars)
	B = numvars["B"]
	A = numvars["A"]
	t1 = 4.0 * pi * G  * A /B *	(2.0/(B**2 * r) - exp(-(r*B)) *  (2.0/(r * B**2) + 2.0 / B  + r))
	pr1 = t1+ K[0] 
	if(pr1<0):
		print("in calc_exp: vc**2 negative=%4.2f t1=%4.2f, K=%4.2f, return 0"%(pr1, t1,K[0]))
		return 0
	return  sqrt(pr1/r)	


def calculateM(numvars, K , r):
	testNumvars(numvars)
	B = numvars["B"]
	A = numvars["A"]
	t1 = e ** (-B * r)
	return 4 * pi * A / B * ( 2* B**(-2) - 2 * B ** (-2) * t1 - (2.0 / B) * r * t1 - r ** 2 * t1 ) + K[0]


