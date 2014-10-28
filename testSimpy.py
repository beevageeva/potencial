from sympy import integrate, Symbol, exp, diff


r = Symbol('r')
r0 = Symbol('r0')
eps =  Symbol('eps')
x =  Symbol('x')
a =  Symbol('a')
#print(integrate(r**2 * exp(-r/r0), r))
#print("--------------Potencial-----")
##print(integrate(1/r**2 * integrate(r**2 * exp(-r/r0), r), r))
#print(integrate(1/x**2 * integrate(a**2 * exp(-a/r0), (a,0,x)), (x,eps,r)))


print("------------vc--")
#print(diff(r0**2*exp(-r/r0) + 2*r0**3*exp(-r/r0)/r, r))
#print(diff(r*integrate(1/x**2 * integrate(a**2 * exp(-a/r0), (a,0,x)), (x,eps,r))  , r))
print(1/r * integrate(a**2 * exp(-a/r0), (a,0,r)))

