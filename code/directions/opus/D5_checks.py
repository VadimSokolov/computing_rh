# D5 experiment 3: tests of the proved bounds and of the failure of the sum reduction. mpmath/numpy floating point, not certified.
import numpy as np, math, mpmath as mp
rng=np.random.default_rng(7)
def mgf_coeffs(w,L):  # symmetric law on -L..L from half weights
    return np.concatenate([w[:0:-1]/2,[w[0]],w[1:]/2])
def zeros_h(p):
    r=np.roots(p[::-1]); return np.log(r.astype(complex))
def strip(p): return np.abs(zeros_h(p).real).max() if len(p)>1 else 0.0
def Y_pair(p1,p2,J):
    L1=(len(p1)-1)//2; L2=(len(p2)-1)//2
    c=np.zeros(2*L1+2*L2+1)
    for i,a in enumerate(p1):
        for j,b in enumerate(p2): c[i+j]+=a*b*math.exp(J*(i-L1)*(j-L2))
    c/=np.abs(c).max(); return np.abs(np.log(np.abs(np.roots(c[::-1])))).max()
print("(a) lattice partner theorem: Y <= max(d2, sqrt(max(d1^2-J^2 M kappa,0))), kappa=1/2 on integer lattice, M=L2 (top mass>0)")
worst=-9; cnt=0
for it in range(4000):
    L1,L2=rng.integers(1,4),rng.integers(1,4)
    p1=mgf_coeffs(rng.dirichlet(np.ones(L1+1)*rng.choice([.2,1])),L1); p2=mgf_coeffs(rng.dirichlet(np.ones(L2+1)*rng.choice([.2,1])),L2)
    d1,d2=strip(p1),strip(p2); J=rng.choice([.02,.1,.3,.7,1.5])
    bd=max(d2,math.sqrt(max(d1*d1-J*J*L2*0.5,0))); y=Y_pair(p1,p2,J)
    worst=max(worst,y-bd); cnt+=1
print(f"  tests={cnt} worst Y-bound={worst:.3e}")
print("(b) identical, split bound: Y <= max(tau, sqrt(d^2-J^2 v_tau)), v_tau=n_tau/4, n_tau=#zeros per period with |Re|<=tau")
worst=-9; cnt=0; gain=0
for it in range(4000):
    L=rng.integers(2,5); p=mgf_coeffs(rng.dirichlet(np.ones(L+1)*rng.choice([.2,.5,1])),L); z=zeros_h(p); re=np.sort(np.abs(z.real)); d=re[-1]
    if d<1e-3: continue
    J=rng.choice([.02,.1,.3,.7]); y=Y_pair(p,p,J); best=d
    for tau in re[:-1]:
        if tau<d-1e-9:
            nt=(np.abs(z.real)<=tau+1e-12).sum(); best=min(best,max(tau,math.sqrt(max(d*d-J*J*nt/4,0))))
    worst=max(worst,y-best); cnt+=1; gain+= best<d-1e-9
print(f"  tests={cnt} worst Y-bound={worst:.3e}, bound strictly below d in {gain} tests")
print("(c) sum reduction: strip of G_J for a zero of multiplicity m (delta=0), y^2/J from E(u+eta^2)^m, eta~N(0,J/2); de Bruijn step removes 1/2")
mp.mp.dps=60
for m in (1,2,3,4,5,6,8,10,16,32):
    co=[mp.binomial(m,k)*mp.fac2(2*k-1)/mp.mpf(2)**k if k>0 else mp.mpf(1) for k in range(m+1)]
    rt=mp.polyroots(co,maxsteps=2000,extraprec=400); v=max((mp.re(u)+abs(u))/2 for u in rt)
    print(f"  m={m:3d} max y^2/J={mp.nstr(v,6)}")
print("(c') law: sigma=sum of 2 fair signs (delta=0), G_J has cosh 2h=-exp(-J) +- i(1-exp(-2J))/sqrt2; y^2/J:")
for J in (1e-3,1e-2,0.1,0.5):
    C=complex(-math.exp(-J),(1-math.exp(-2*J))/math.sqrt(2)); h=np.arccosh(C+0j)/2
    print(f"  J={J}: Re h={abs(h.real):.6f}, (Re h)^2/J={h.real**2/J:.4f}")
print("(d) sharpness: Ising partner lambda=1, first spin psi^m, psi=spin one scale omega with strip delta; zero near h0=y*+i pi/omega")
mp.mp.dps=40
for (omega,delta,J) in ((1.0,0.9624,0.3),(1.0,0.9624,0.7),(0.1,1.0,0.5),(0.01,1.0,0.5)):
    ys=mp.acosh(mp.cosh(omega*delta)/mp.cosh(omega*J))/omega
    psi=lambda h:mp.cosh(omega*delta)+mp.cosh(omega*h)
    row=[]
    for m in (3,11,101,1001,10001):
        F=lambda h: psi(h+J)/psi(h-J)+mp.exp(-2*h/m)
        h=mp.findroot(F,mp.mpc(ys,mp.pi/omega))
        row.append(mp.re(h))
    print(f"  omega={omega} delta={delta} J={J}: y*={mp.nstr(ys,6)}, sqrt(d^2-J^2)={math.sqrt(delta**2-J*J):.6f}, Re h_m(m=3,11,101,1001,10001)="+" ".join(mp.nstr(x,6) for x in row))
