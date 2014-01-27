from scipy import integrate
from math import e, sqrt, pi
from numpy import inf
import sys
#for the potential there will be two integration constants K[0] and K[1]
#for the rest(velocity , masss, integrated mass) there is only one K[0]
#numvars is the hash containing variable parameters(A, B, ...)

G = 6.6 * 10 **(-11)

def testNumvars(numvars):
	if(not 'A' in numvars.keys()  ):
		print("parameter A  not present in numvars:")
		print(numvars)
		sys.exit(0)
	

def calculateP(numvars,K,r, startPoint):
	testNumvars(numvars)
	A = numvars["A"]
	#print("A=%4.2f, startpoint=%4.2f"%(A,startPoint))
	int1 = integrate.quad(lambda a: a ** (-2) * e ** (-A * a) , startPoint, r )[0]
	int2 = integrate.quad(lambda a: a ** (-1) * e ** (-A * a) , startPoint, r )[0]
	result = (4 * pi * G  / A)* (-(2.0  * int1 / A) - 2.0  * int2  +  e ** (-A * r)) + K[0] / r  + K[1]
	return result

def calculateV(numvars, K , r, startPoint):
	testNumvars(numvars)
	A = numvars["A"]
	t1 =  (-4.0 * pi * G ) * (2 * A **(-2) * r **(-2) * e ** (-A * r) + 2 * A ** (-1) * r ** (-1) * e ** (-A * r) +  e **(-A * r) )	
	pr1 = t1+ K[0] 
	if(pr1<0):
		print("in calc_exp: vc**2 negative=%4.2f t1=%4.2f, K=%4.2f, return 0"%(pr1, t1,K[0]))
		return 0
	return  sqrt(pr1/r)	


def calculateM(numvars, K , r, startPoint):
	testNumvars(numvars)
	A = numvars["A"]
	return (-4 * pi ) * (2 * A ** (-2) * e ** (-A * r) + 2 / A * r * e ** (-A * r) + r ** 2 * e ** (-A * r) ) + K[0]

def calculateDp(numvars, K , s, startPoint):
	testNumvars(numvars)
	A = numvars["A"]
	int1 =  integrate.quad(lambda r: (r * e ** (-A * r)) / sqrt(r**2-s**2)  , s, inf)
	return 2 * int1[0]	

