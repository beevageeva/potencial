from sympy import Symbol, solve


a = Symbol("a")
b = Symbol("b")
print(solve(a-b, a))
