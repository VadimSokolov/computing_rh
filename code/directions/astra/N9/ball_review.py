"""Fixed absolute precision dyadic midpoint-radius balls for this review.

Arithmetic uses exact Python integers with outward rounding. Elementary
functions and constants are enclosed by mpmath.iv, then rounded outwards.
This is an independent small certificate implementation, not Arb.
"""
from fractions import Fraction
import mpmath as mp

BITS = 256
SCALE = 1 << BITS
mp.iv.prec = BITS + 32


def set_precision(bits):
    global BITS, SCALE
    BITS = bits
    SCALE = 1 << bits
    mp.iv.prec = bits + 32


def _floor_endpoint(t):
    sign, man, exp, bc = t
    if bc < 0:
        raise ArithmeticError('nonfinite interval endpoint')
    man = -man if sign else man
    exp += BITS
    return man << exp if exp >= 0 else man >> (-exp)


def _ceil_endpoint(t):
    sign, man, exp, bc = t
    return -_floor_endpoint((not sign, man, exp, bc))


class B:
    __slots__ = ('m', 'r')

    def __init__(self, value=0):
        if isinstance(value, B):
            self.m, self.r = value.m, value.r
            return
        f = Fraction(value)
        num = f.numerator*SCALE
        self.m = num//f.denominator
        self.r = 0 if num % f.denominator == 0 else 1

    @classmethod
    def raw(cls, m, r):
        obj = object.__new__(cls)
        obj.m, obj.r = int(m), int(r)
        return obj

    @classmethod
    def from_iv(cls, x):
        lo = _floor_endpoint(x._mpi_[0])
        hi = _ceil_endpoint(x._mpi_[1])
        mid = (lo+hi)//2
        return cls.raw(mid, max(mid-lo,hi-mid))

    def interval(self):
        return mp.iv.mpf([mp.iv.mpf(self.m-self.r)/SCALE, mp.iv.mpf(self.m+self.r)/SCALE])

    def __add__(self, other):
        other=B(other)
        return B.raw(self.m+other.m,self.r+other.r)

    __radd__=__add__

    def __neg__(self):
        return B.raw(-self.m,self.r)

    def __sub__(self, other):
        return self+-B(other)

    def __rsub__(self, other):
        return B(other)+-self

    def __mul__(self, other):
        other=B(other)
        product=self.m*other.m
        mid=product//SCALE
        err=abs(self.m)*other.r+abs(other.m)*self.r+self.r*other.r
        rad=(err+SCALE-1)//SCALE+(product % SCALE != 0)
        return B.raw(mid,rad)

    __rmul__=__mul__

    def __truediv__(self, other):
        other=B(other)
        low,high=other.m-other.r,other.m+other.r
        if low<=0<=high:
            raise ArithmeticError('division by ball containing zero')
        vals=[(x*SCALE,y) for x in [self.m-self.r,self.m+self.r] for y in [low,high]]
        lo=min(x//y for x,y in vals)
        hi=max(-((-x)//y) for x,y in vals)
        mid=(lo+hi)//2
        return B.raw(mid,max(mid-lo,hi-mid))

    def __rtruediv__(self, other):
        return B(other)/self

    def positive(self):
        return self.m-self.r>0

    def negative(self):
        return self.m+self.r<0

    def bounds(self):
        return ((self.m-self.r)/SCALE,(self.m+self.r)/SCALE)

    def mid(self):
        return self.m/SCALE

    def radius(self):
        return self.r/SCALE

    def __repr__(self):
        return f'{self.mid():.15g} +/- {self.radius():.3g}'


def blog(x):
    return B.from_iv(mp.iv.log(B(x).interval()))


def bexp(x):
    return B.from_iv(mp.iv.exp(B(x).interval()))


def bsqrt(x):
    return B.from_iv(mp.iv.sqrt(B(x).interval()))


def ball_levinson(c):
    coeff=[]
    v=c[0]
    vs=[v]
    ks=[]
    if not v.positive():
        return vs,ks,1
    for k in range(1,len(c)):
        dot=sum((coeff[j]*c[k-j-1] for j in range(k-1)),B(0))
        refl=(c[k]-dot)/v
        ks.append(refl)
        coeff=[coeff[j]-refl*coeff[-j-1] for j in range(k-1)]+[refl]
        v=v*(1-refl*refl)
        vs.append(v)
        if not v.positive():
            return vs,ks,k+1
    return vs,ks,None
