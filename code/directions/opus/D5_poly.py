# D5 experiment 1: identical spins, local polynomial model.
# Q(w) has zeros symmetric under r -> -conj(r), all in |Re r| <= 1.  Z(w) = sum_k J^k/k! (Q^(k)(w))^2.
# Compares ymax^2 with 1 - J/2 and with f(J) = (1 - 2J + sqrt(1+2J^2))/2 (exact for one pair).
# Floating point (numpy.roots), not certified.
import numpy as np, math, sys
P=np.polynomial.polynomial
def Zpoly(q,J):            # q: ascending coefficients
    out=np.zeros(2*len(q)-1,dtype=complex); d=q.copy(); fac=1.0; k=0
    while True:
        m=P.polymul(d,d); out[:len(m)]+=fac*m
        if len(d)==1: break
        d=P.polyder(d); k+=1; fac*=J/k
    return out
def ymax(q,J):
    r=P.polyroots(Zpoly(q,J)); return r.real.max(),r
def Qfrom(ts,ds):
    q=np.array([1.0+0j])
    for t,d in zip(ts,ds):
        q=P.polymul(q,np.array([(1j*t)**2-d**2,2j*t,1.0]))   # (w+it)^2-d^2
    return q
f=lambda J:(1-2*J+math.sqrt(1+2*J*J))/2 if J<2 else 0.0
rng=np.random.default_rng(5)
Js=[0.02,0.05,0.1,0.2,0.3,0.5,0.75,1.0,1.25,1.5,1.75,1.9,2.0,2.2]
print("J      f(J)    1-J/2  | n: max ymax^2 over random heights (all pairs at Re=+-1)")
for J in Js:
    row=[]
    for n in (1,2,3,4):
        best=-1
        for it in range(1 if n==1 else 4000):
            s=rng.choice([0.1,0.3,1,3,10])
            ts=rng.normal(0,s*math.sqrt(max(J,0.01)),n); ts[0]=0
            y,_=ymax(Qfrom(ts,np.ones(n)),J)
            if y>best: best=y
        row.append(max(best,0)**2)
    print(f"{J:5.2f} {f(J):7.4f} {max(1-J/2,0):7.4f} | "+" ".join(f"{x:7.4f}" for x in row))
