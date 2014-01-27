from scipy import integrate
from math import e, sqrt, pi
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
	

def calculateP(numvars,K,r, startPoint):
	testNumvars(numvars)
	B = numvars["B"]
	A = numvars["A"]
	int1 = integrate.quad(lambda a: a ** (-2) * e ** (-B * a) , startPoint, r )[0]
	int2 = integrate.quad(lambda a: a ** (-1) * e ** (-B * a) , startPoint, r )[0]
	result = (4 * pi * G * A / B ** 2)* (-(2.0  * int1 / B) - 2.0  * int2  +  e ** (-B * r)) + K[0] / r  + K[1]
	return result

def calculateV(numvars, K , r, startPoint):
	testNumvars(numvars)
	B = numvars["B"]
	A = numvars["A"]
	t1 =  (-4.0 * pi * G  * A *  B ** (-1)) * (2 * B **(-2) * r **(-2) * e ** (-B * r) + 2 * B ** (-1) * r ** (-1) * e ** (-B * r) +  e **(-B * r) )	
	pr1 = t1+ K[0] 
	if(pr1<0):
		print("in calc_exp: vc**2 negative=%4.2f t1=%4.2f, K=%4.2f, return 0"%(pr1, t1,K[0]))
		return 0
	return  sqrt(pr1/r)	


def calculateM(numvars, K , r, startPoint):
	testNumvars(numvars)
	B = numvars["B"]
	A = numvars["A"]
	return (-4 * pi * A / B) * (2 * B ** (-2) * e ** (-B * r) + 2 / B * r * e ** (-B * r) + r ** 2 * e ** (-B * r) ) + K[0]

def calculateDp(numvars, K , s, startPoint):
	testNumvars(numvars)
	B = numvars["B"]
	A = numvars["A"]
	int1 =  integrate.quad(lambda r: (r * e ** (-B * r)) / sqrt(r**2-s**2)  , s, inf)
	return 2 * int1[0]	

