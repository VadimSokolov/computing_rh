# D5 experiment 4: two Polya spins with a planted pair (k=1, eta=0.3, as in the book), pure pair coupling J versus Curie Weiss a.
# Zeros r of L(h)=xi(1/2+h)/xi(1/2): +-i gamma_k (k>=3) and +-eta +- i(gamma1+gamma2)/2; truncated at K with Gaussian tail exp(B h^2).
# Z(h,h)=sum_k J^k/k! (L^(k)(h))^2 via Taylor coefficients at h. Double precision, not certified. The planted function is not known to be a law.
import numpy as np, math, mpmath as mp
K=150; N=36; eta=0.3
gam=np.array([float(mp.zetazero(n).imag) for n in range(1,K+1)])
xi=lambda s:0.5*s*(s-1)*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
V=float(mp.diff(lambda s:mp.log(xi(s)),mp.mpf(1)/2,2))      # =2*sum 1/gamma^2
B=V/2-(1/gam**2).sum()
g0=(gam[0]+gam[1])/2; r=complex(eta,g0)
fac=np.array([math.factorial(k) for k in range(N+1)],dtype=float)
def far(planted):
    if planted:
        z=np.concatenate([1j*gam[2:],-1j*gam[2:],[complex(eta,-g0),complex(-eta,-g0)]]); return z,np.array([complex(eta,g0),complex(-eta,g0)])
    return np.concatenate([1j*gam[1:],-1j*gam]),np.array([1j*gam[0]])
def coeffs(h,planted):
    zf,zn=far(planted); d=h-zf
    a=np.zeros(N+1,dtype=complex); pw=np.ones_like(d)
    for n in range(1,N+1):
        pw=pw/d; a[n]=(-1)**(n-1)*pw.sum()/n
    a[1]+=2*B*h; a[2]+=B
    c=np.zeros(N+1,dtype=complex); c[0]=1
    for k in range(1,N+1):
        c[k]=sum(j*a[j]*c[k-j] for j in range(1,k+1))/k
    for z in zn:
        c=np.convolve(c,[h-z,1.0])[:N+1]
    return c
def Zpair(h,J,planted=True):
    c=coeffs(h,planted); return np.sum(fac*J**np.arange(N+1)*c*c)
def Zcw(h,a,planted=True):
    c=coeffs(h,planted); d=np.convolve(c,c)[:N+1]; s=0j; w=1.0
    for k in range(0,N//2+1):
        s+=w*d[2*k]; w*=a*(2*k+1)
    return s
def secant(F,x0):
    x1=x0+1e-4*(1+1j); f0,f1=F(x0),F(x1)
    for _ in range(60):
        if f1==f0: break
        x2=x1-f1*(x1-x0)/(f1-f0); x0,f0=x1,f1; x1=x2; f1=F(x1)
        if abs(x1-x0)<1e-12: break
    return x1
# first order drift of the rightmost zero: sum_{r'!=r} 1/(r-r')
zf,zn=far(True); c1=(1/(r-zf)).sum()+1/(r-zn[1])+2*B*r
print(f"V={V:.6f} B(tail)={B:.3e} r={r:.5f}; Re sum 1/(r-r')={c1.real:.4f} (mirror alone 1/(2 eta)={1/(2*eta):.4f}); first order rate d(delta^2)/dJ=-{2*eta*c1.real:.3f}")
print("Curie Weiss check against the book (0.293,0.265,0.229,0.149 at a=.002,.01,.02,.04):")
for a in (0.002,0.01,0.02,0.04):
    u=secant(lambda h:Zcw(h,a),r+1j*math.sqrt(2*a)-2*a*c1); print(f"  a={a}: Re h={u.real:.4f} Im h-g0={u.imag-g0:.4f}")
f=lambda u:(1-2*u+math.sqrt(1+2*u*u))/2 if u<2 else 0.0
print("pure pair coupling J: branches from r +- i sqrt(J)")
print(" J       Re(up)   Re(low)  (eta^2-max Re^2)/J   sqrt(eta^2-J)    one pair sqrt(eta^2 f)")
up,lo=r+0j,r+0j
for J in (0.001,0.002,0.005,0.01,0.02,0.03,0.04,0.06,0.08,0.1,0.14,0.18,0.25):
    if J<=0.002: up,lo=r+1j*math.sqrt(J)-J*c1,r-1j*math.sqrt(J)-J*c1
    up=secant(lambda h:Zpair(h,J),up); lo=secant(lambda h:Zpair(h,J),lo)
    m=max(abs(up.real),abs(lo.real))
    print(f"{J:6.3f} {up.real:8.4f} {lo.real:8.4f} {(eta*eta-m*m)/J:10.3f} {math.sqrt(max(eta*eta-J,0)):14.4f} {eta*math.sqrt(f(J/eta**2)):14.4f}")
print("true Polya spins, branches from i gamma1 +- i sqrt(J):")
for J in (0.001,0.01,0.05):
    a=secant(lambda h:Zpair(h,J,False),1j*gam[0]+1j*math.sqrt(J)); b=secant(lambda h:Zpair(h,J,False),1j*gam[0]-1j*math.sqrt(J))
    print(f"  J={J}: Re={a.real:.2e},{b.real:.2e} Im-gamma1={a.imag-gam[0]:.5f},{b.imag-gam[0]:.5f}")
