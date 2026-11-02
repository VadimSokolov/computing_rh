# Davenport-Heilbronn: u_alpha^DH(w)=(1/pi) Re Lambda'/Lambda(alpha+iw). Own code.
import numpy as np
from scipy.special import digamma
from math import factorial
B=[1/6,-1/30,1/42,-1/30,5/66,-691/2730,7/6]
kap=(np.sqrt(10-2*np.sqrt(5))-2)/(np.sqrt(5)-1)
coef=[(1,1.0),(2,kap),(3,-kap),(4,-1.0)]
def hz(s,a,N=400):
    n=np.arange(N)+a; ln=np.log(n); t=np.exp(-s*ln); z=t.sum(); dz=-(ln*t).sum()
    M=N+a; LN=np.log(M); P=np.exp(-s*LN)
    z+=P*M/(s-1)+.5*P; dz+=-LN*P*M/(s-1)-P*M/(s-1)**2-.5*LN*P
    poch=s;dp=1+0j;term=P/M
    for k,b in enumerate(B,1):
        c=b/factorial(2*k)
        z+=c*poch*term; dz+=c*(dp*term-LN*poch*term)
        a1=s+2*k-1;a2=s+2*k
        dp=dp*a1*a2+poch*(a1+a2);poch=poch*a1*a2;term=term/M/M
    return z,dz
def f(s):
    z=0;dz=0
    for k,c in coef:
        a,d=hz(s,k/5); z+=c*a; dz+=c*d
    q=5**(-s)
    return q*z, q*dz-np.log(5)*q*z
def dL(s):
    z,dz=f(s)
    return .5*np.log(5/np.pi)+.5*digamma((s+1)/2)+dz/z
print("Re at 1/2",[f"{dL(.5+1j*w).real:.1e}" for w in (3.,20.,85.)])
print("|f(.808517+85.699348i)|",abs(f(.808517+85.699348j)[0]))
w0=85.699348
for al in [.51,.55,.6,.65,.7,.75,.78,.8,.805,.82,.9,1.0,1.1,1.2,1.3,1.5]:
    W=np.arange(0,120.0001,.05)
    u=np.array([dL(al+1j*w).real/np.pi for w in W])
    print(f"alpha {al}: min u on [0,120] {u.min():+.4f} at w={W[u.argmin()]:.2f}; u(85.699)={dL(al+1j*w0).real/np.pi:+.4f}; neg on [80,92] frac {np.mean(u[(W>80)&(W<92)]<0):.3f}; neg overall {np.mean(u<0):.3f}")
# large scan for alpha>1: is there negativity? zero beyond 1 appears at much larger height
for al in [1.0,1.05,1.1,1.2]:
    W=np.arange(0,4000,.1)
    u=[dL(al+1j*w).real/np.pi for w in W[::1]]
    u=np.array(u)
    print(f"alpha {al}: min on [0,4000] {u.min():+.4f} at {W[u.argmin()]:.1f}")
# kernel positivity at x>=1
x=np.linspace(1,20,2000)
g=sum(c*n*np.exp(-np.pi*n*n*x/5) for n,c in [(1,1),(2,kap),(3,-kap),(4,-1),(6,1),(7,kap),(8,-kap),(9,-1)])
print("kernel min x>=1 scaled",(g*np.exp(np.pi*x/5)).min())
