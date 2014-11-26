from sympy import integrate, Symbol, exp, diff, sqrt


r = Symbol('r')
b = Symbol('b')
M =  Symbol('M')
G =  Symbol('G')
pi =  Symbol('pi')

#isochrone
pot = - (G * M) / (b + sqrt(b**2 + r**2))
#b = 0, all mass in 0
#pot = - (G * M) / r



dif1Pot =  diff(pot  , r)

print("------------vc-------------")
print(sqrt(r * dif1Pot))

print("------------mass-------------")
print(1/G * r**2 * dif1Pot)


print("------------dens-------------")
print(1 / (4 * pi * G) * 1 / (r**2) * diff(r**2 * dif1Pot, r))

