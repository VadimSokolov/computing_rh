# Model test: xi with two adjacent on-line zeros replaced by an off-line quadruple.
# Euler-Maclaurin in double precision, N >= T. Own code.
import numpy as np, sys
from scipy.special import digamma, loggamma
B=[1/6,-1/30,1/42,-1/30,5/66,-691/2730,7/6]
def zz(s,N):
    n=np.arange(1,N,dtype=float); ln=np.log(n)
    t=np.exp(-s*ln); z=t.sum(); dz=-(ln*t).sum()
    LN=np.log(N); P=np.exp(-s*LN)
    z+=P*N/(s-1)+.5*P
    dz+=-LN*P*N/(s-1)-P*N/(s-1)**2-.5*LN*P
    # Bernoulli tail, value and s-derivative
    poch=s; dpoch=1.0+0j; fact=2.0; Np=N**-1.0
    term=P/N
    from math import factorial
    for k,b in enumerate(B,1):
        c=b/factorial(2*k)
        z+=c*poch*term
        dz+=c*(dpoch*term-LN*poch*term)
        # next (s)_{2k+1}
        a1=s+2*k-1;a2=s+2*k
        dpoch=dpoch*a1*a2+poch*(a1+a2); poch=poch*a1*a2; term=term/N/N
    return z,dz
def dlxi(s,N):
    z,dz=zz(s,N)
    return 1/s+1/(s-1)-.5*np.log(np.pi)+.5*digamma(s/2)+dz/z
def Z(w,N):
    s=.5+1j*w; z,_=zz(s,N)
    th=loggamma(.25+.5j*w).imag-.5*w*np.log(np.pi)
    return (np.exp(1j*th)*z).real
assert abs(zz(2+0j,2000)[0]-np.pi**2/6)<1e-10
for T in [1e4,1e5,1e6]:
    N=int(max(T,3000))
    # locate zeros near T
    ws=np.arange(T,T+14,.05); zv=[Z(w,N) for w in ws]
    zs=[]
    for i in range(len(ws)-1):
        if zv[i]*zv[i+1]<0:
            a,b=ws[i],ws[i+1];fa=zv[i]
            for _ in range(30):
                m=.5*(a+b);fm=Z(m,N)
                if fa*fm<=0:b=m
                else:a,fa=m,fm
            zs.append(.5*(a+b))
    zs=np.array(zs)
    i=len(zs)//2; g1,g2=zs[i],zs[i+1]; gb=.5*(g1+g2)
    print("T",T,"zeros in window",len(zs),"pair",round(g1,5),round(g2,5))
    for al,bet in [(.75,.95),(.75,.9),(.75,.85),(.75,.8),(.6,.95),(.6,.8),(.9,.98)]:
        def uG(w):
            s=al+1j*w
            v=dlxi(s,N)
            rho=bet+1j*gb
            for r in [rho,bet-1j*gb,1-bet+1j*gb,1-bet-1j*gb]: v+=1/(s-r)
            for g in [g1,g2,-g1,-g2]: v-=1/(s-.5-1j*g)
            return v.real/np.pi
        grid=np.arange(gb-6,gb+6.0001,.05)
        vals=np.array([uG(w) for w in grid])
        print(f"  alpha {al} beta {bet} eta {bet-al:.2f}: min u^G on +-6 = {vals.min():+.4f} at w-gb={grid[vals.argmin()]-gb:+.2f}; u^G(gb)={uG(gb):+.4f}")
