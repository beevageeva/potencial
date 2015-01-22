from scipy import integrate as integrate_num
import numpy as np
from math import pi
from sympy import oo,Symbol, integrate

int1 = integrate_num.quad(lambda x: (x**2+1)**(-3/2), -np.inf, np.inf)
print("NUMERIC")
print(int1[0])
print("pi/2=%1.5f" % (pi/2.0))
#analytical
print("ANALITIC")
x = Symbol('x')
print(integrate((x**2+1)**(-3/2),x))
print(integrate((x**2+1)**(-3/2),(x, -oo,oo)))
