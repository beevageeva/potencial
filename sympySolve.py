from sympy import exp, sqrt, Symbol, solve

M = Symbol("M")
pi = Symbol("pi")
G = Symbol("G")
Rmax = Symbol("Rmax")
r0 = Symbol("r0")
Rsun = Symbol("Rsun")
vSun = Symbol("vSun")
rhoC = Symbol("rhoC")

#not implemented error
#t1 = exp(-(Rmax/r0))
t1 = 1 - Rmax/r0 + 1/2 * Rmax**2/r0**2 - 1/6 * Rmax**3/r0**3

#print(solve([4 * pi * rhoC * r0 * (2* r0**2 - 2 * r0 **2 * t1 - 2.0 *r0 * Rmax * t1 - Rmax ** 2 * t1 )-M ,  sqrt(4.0 * pi * G  * rhoC *r0* (2 *  r0**2 / Rsun - exp(-(Rsun/r0)) *  (2 * r0**2/Rsun + 2 * r0  + Rsun))) - vSun], r0, rhoC, check=False, numerical=False,  dict=True, set=True))

#not implemented error
print(solve([4 * pi * rhoC * r0 * (2* r0**2 - 2 * r0 **2 * t1 - 2.0 *r0 * Rmax * t1 - Rmax ** 2 * t1 )-M ,  sqrt(4.0 * pi * G  * rhoC *r0* (2 *  r0**2 / Rsun - (1-Rsun/r0 + 1/2 * Rsun**2/r0**2 - 1/6 * Rsun**3/r0**3 ) *  (2 * r0**2/Rsun + 2 * r0  + Rsun))) - vSun], r0, rhoC, check=False, numerical=False,  dict=True, set=True))
