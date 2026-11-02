# D5 experiment 7: open chains of N identical spins, local pair Q(w)=w^2-1, common field, exp(J sum d_i d_{i+1}). Numpy roots, not certified.
# Perturbation prediction: d(y^2)/dJ at 0 equals -min over path modes of <d> = -(2-(4/(N+1)) max_k sin^2(k pi/(N+1))) (see D5.md).
import numpy as np, math, itertools
P=np.polynomial.polynomial
Qd=[np.array([-1.0,0,1.0]),np.array([0,2.0]),np.array([2.0])]   # Q, Q', Q'' ascending
def chainZ(N,J):
    out=np.zeros(2*N+1)
    for ks in itertools.product((0,1,2),repeat=N-1):
        ord=[0]*N
        for e,k in enumerate(ks): ord[e]+=k; ord[e+1]+=k
        if max(ord)>2: continue
        coef=1.0
        for k in ks: coef*=J**k/math.factorial(k)
        p=np.array([1.0])
        for o in ord: p=P.polymul(p,Qd[o])
        out[:len(p)]+=coef*p
    return out
def ymax(N,J):
    c=chainZ(N,J); c=np.trim_zeros(c,'b'); r=P.polyroots(c); return max(r.real.max(),0.0)
def rate(N):
    return 2-(4/(N+1))*max(math.sin(k*math.pi/(N+1))**2 for k in range(1,N+1))
print("N  predicted rate | measured (1-y^2)/J at J=1e-4,1e-3 | y^2 at J=.1,.25,.5,1,1.5 | closure J*")
for N in (2,3,4,5,6,8):
    r1=(1-ymax(N,1e-4)**2)/1e-4; r2=(1-ymax(N,1e-3)**2)/1e-3
    row=[ymax(N,J)**2 for J in (0.1,0.25,0.5,1,1.5)]
    lo,hi=0.0,4.0
    for _ in range(30):
        mid=(lo+hi)/2
        if ymax(N,mid)>1e-6: lo=mid
        else: hi=mid
    print(f"{N:2d} {rate(N):8.4f} | {r1:8.4f} {r2:8.4f} | "+" ".join(f"{x:.4f}" for x in row)+f" | {hi:.4f}")
