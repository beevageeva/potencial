from sympy import integrate, Symbol, exp, sqrt


r = Symbol('r')
eps =  Symbol('eps')
x =  Symbol('x')
a =  Symbol('a')
pi = Symbol('pi')
G = Symbol('G')

#define density function
r0 = Symbol('r0')
rhoC = Symbol('rhoC')
densA = rhoC * exp(-a/r0)


int1X = integrate(a**2 * densA, (a,0,x))

print("--------------Potencial(r)-----")
print(4 * pi * G * integrate(1/x**2 * int1X, (x,eps,r)))


print("------------vc(x)----------------")
print(sqrt((4 * pi * G)/x *int1X ))

print("------------mass(x)----------------")
print(4 * pi *int1X )
