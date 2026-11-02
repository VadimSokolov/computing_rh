from decimal import Decimal, Context, ROUND_FLOOR, ROUND_CEILING, ROUND_HALF_EVEN
from fractions import Fraction
from math import factorial

PREC = 140
dn = Context(prec=PREC, rounding=ROUND_FLOOR, Emin=-9999999, Emax=9999999)
up = Context(prec=PREC, rounding=ROUND_CEILING, Emin=-9999999, Emax=9999999)
rn = Context(prec=PREC, rounding=ROUND_HALF_EVEN, Emin=-9999999, Emax=9999999)

class I:
    def __init__(self, lo, hi=None):
        self.lo = Decimal(lo)
        self.hi = self.lo if hi is None else Decimal(hi)
        assert self.lo <= self.hi
    @staticmethod
    def cast(x):
        return x if isinstance(x, I) else I(x)
    def __add__(self, other):
        other = I.cast(other)
        return I(dn.add(self.lo, other.lo), up.add(self.hi, other.hi))
    __radd__ = __add__
    def __neg__(self):
        return I(self.hi.copy_negate(), self.lo.copy_negate())
    def __sub__(self, other):
        return self + (-I.cast(other))
    def __rsub__(self, other):
        return I.cast(other) - self
    def __mul__(self, other):
        other = I.cast(other)
        pairs = [(x,y) for x in (self.lo,self.hi) for y in (other.lo,other.hi)]
        return I(min(dn.multiply(x,y) for x,y in pairs), max(up.multiply(x,y) for x,y in pairs))
    __rmul__ = __mul__
    def __truediv__(self, other):
        other = I.cast(other)
        assert not (other.lo <= 0 <= other.hi)
        reciprocal = I(dn.divide(Decimal(1),other.hi),up.divide(Decimal(1),other.lo))
        return self*reciprocal
    def __pow__(self, n):
        assert isinstance(n,int) and n >= 0
        answer, base = I(1), self
        while n:
            if n & 1:
                answer = answer*base
            base = base*base
            n //= 2
        return answer
    def exp(self):
        # Decimal.exp is correctly rounded to nearest. One adjacent value
        # on each side therefore supplies an outward enclosure.
        return I(rn.next_minus(rn.exp(self.lo)),rn.next_plus(rn.exp(self.hi)))
    def show(self):
        return '['+str(self.lo)+', '+str(self.hi)+']'

def rational(x):
    x = Fraction(x)
    return I(dn.divide(Decimal(x.numerator),Decimal(x.denominator)),up.divide(Decimal(x.numerator),Decimal(x.denominator)))

def atan_inverse(n):
    terms = 128
    s = sum((Fraction((-1)**k,(2*k+1)*n**(2*k+1)) for k in range(terms)), Fraction(0))
    remainder = Fraction(1,(2*terms+1)*n**(2*terms+1))
    return rational(s)+I(0,rational(remainder).hi)

pi = 16*atan_inverse(5)-4*atan_inverse(239)
h = rational(Fraction(1,128))
degree = 14
raw = [I(0) for _ in range(degree+1)]
for k in range(513):
    x = k*h
    t = (2*x).exp()
    t2 = t*t
    kernel = I(0)
    for n in range(1,13):
        kernel += (4*pi*pi*n**4*t2-6*pi*n*n*t)*(-pi*n*n*t).exp()
    kernel *= (x/2).exp()
    weight = h*(1 if k == 0 else 2)*kernel
    power = I(1)
    for m in range(degree+1):
        raw[m] += weight*power
        power *= x*x
error = I('-1e-120','1e-120')
raw = [v+error for v in raw]
mom = [v/raw[0] for v in raw]
mom[0] = I(1)
b = [mom[m]/factorial(2*m) for m in range(degree+1)]
logb = [I(0) for _ in range(degree+1)]

def bernoulli(n):
    a = [Fraction(0) for _ in range(n+1)]
    for m in range(n+1):
        a[m] = Fraction(1,m+1)
        for j in range(m,0,-1):
            a[j-1] = j*(a[j-1]-a[j])
    return a[0]

T = [I(0) for _ in range(degree+1)]
for m in range(1,degree+1):
    logb[m] = b[m]-sum((k*logb[k]*b[m-k] for k in range(1,m)),I(0))/m
    kappa = factorial(2*m)*logb[m]
    uniform = Fraction(2**(2*m),2*m)*bernoulli(2*m)
    T[m] = 20**m*kappa/rational(uniform)

q = [I(s) for s in ['-0.0000840961934787','0.00209340425559','0.00472436230707','-0.222651479574','1.08863834345','-1.86877263496','1']]
answer = sum((q[i]*q[j]*T[i+j+2] for i in range(7) for j in range(7)),I(0))
print('pi',pi.show())
print('mass',raw[0].show())
print('quadratic_form',answer.show())
print('width',up.subtract(answer.hi,answer.lo))
assert answer.hi < Decimal('-6.8159099165e-11')
assert answer.lo > Decimal('-6.8159099166e-11')
for m in range(2,15):
    print('T',m,T[m].show())
print('INTERVAL ENCLOSURE PASSED')
