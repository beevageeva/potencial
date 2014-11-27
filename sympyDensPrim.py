from sympy import integrate, Symbol, exp, sqrt


r = Symbol('r')
eps =  Symbol('eps')
pi = Symbol('pi')
G = Symbol('G')

#define density function
r0 = Symbol('r0')
rhoC = Symbol('rhoC')
densR = rhoC * exp(-r/r0)


int1 = integrate(r**2 * densR,r)

print("--------------Potencial(r)-----")
print(4 * pi * G * integrate(1/r**2 * int1, r))


print("------------vc(x)----------------")
print(sqrt((4 * pi * G)/r *int1 ))

print("------------mass(x)----------------")
print(4 * pi *int1 )
