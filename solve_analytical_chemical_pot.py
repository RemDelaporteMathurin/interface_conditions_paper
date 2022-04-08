from sympy import *

D1, D2 = symbols('D1, D2')
x_int = Symbol("x_int")
S1, S2 = symbols('S1, S2')
u_0, u_L = symbols('c_0, c_L')
a1, a2, b1, b2 = symbols('a1, a2, b1, b2')
x = Symbol("x")
f = Symbol("f")
u1 = 1/2*f/D1*x**2 + a1*x + b1
u2 = 1/2*f/D2*x**2 + a2*x + b2
x_L = Symbol('L')
x_0 = 0#Symbol('x_0')

list_of_equations = [
    u1.subs(x, x_int)/S1 - u2.subs(x, x_int)/S2,
    D1*diff(u1, x) - D2*diff(u2, x),
    u1.subs(x, x_0) - u_0,
    u2.subs(x, x_L) - u_L
]
res = solve(list_of_equations, a1, a2, b1, b2)

print(res[b2])
print(latex(res[b2]))

# print(latex(eval(str(res[b1]))))
# print(latex(eval(str(res[a2]))))
# print(latex(eval(str(res[b2]))))
