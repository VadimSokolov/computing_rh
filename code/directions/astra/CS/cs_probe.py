import sys
sys.path.insert(0, '/tmp/cs_lib')
import flint
from flint import arb, arb_poly, acb, ctx
ctx.prec = 160
print('Python', sys.version.split()[0], 'python-flint', flint.__version__)
t = arb_poly([0, 1])
r = 3 / (2 * arb.pi())
f = 1 - r*t**2 - r*t**7 + t**9
def d(p):
    return t * p.derivative()
h = d(d(d(f)))*f**2 - 3*d(d(f))*d(f)*f + 2*d(f)**3
p = 8*arb.pi()*(t**4-1)*f**3 - t**2*h
q = p(arb_poly([1, 1]))
print('GHS polynomial degree', q.degree())
for j in range(q.degree()+1):
    print(j, q[j])
print('Acb integral doc:', acb.integral.__doc__)
