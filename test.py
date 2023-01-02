from sympy import Symbol

x = Symbol('x')
y = Symbol('y')

z = x * y
w = z.diff(x)
print(w.subs({y:10, x :9}))

