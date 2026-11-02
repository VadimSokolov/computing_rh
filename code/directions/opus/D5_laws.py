# D5 experiment 2: true laws on integer lattices and Ising partners. Floating point, not certified.
import numpy as np, math
rng=np.random.default_rng(55)
P=np.polynomial.polynomial
def strip_law(p):                 # p[s+L], s=-L..L, zeros of mgf in z=e^h
    r=np.roots(p[::-1]); return np.abs(np.log(np.abs(r))).max() if len(r) else 0.0
def Y_pair(p1,p2,J):              # common field, exp(J s1 s2)
    L1=(len(p1)-1)//2; L2=(len(p2)-1)//2
    c=np.zeros(2*L1+2*L2+1)
    for i,a in enumerate(p1):
        for j,b in enumerate(p2):
            c[i+j]+=a*b*math.exp(J*(i-L1)*(j-L2))
    c/=c.max(); r=np.roots(c[::-1]); return np.abs(np.log(np.abs(r))).max()
f=lambda u:(1-2*u+math.sqrt(1+2*u*u))/2 if u<2 else 0.0
print("(1) identical spin one: closure J*=log cosh 2a, ymax at J*-eps and J*+eps")
for q in (0.2,0.1,0.24,0.249):
    p=np.array([q,1-2*q,q]); a=math.acosh((1-2*q)/(2*q)); Js=math.log(math.cosh(2*a))
    print(f" q={q} a={a:.4f} J*={Js:.5f} Y(J*-.001)={Y_pair(p,p,Js-1e-3):.4f} Y(J*+.001)={Y_pair(p,p,Js+1e-3):.2e} 2a^2={2*a*a:.5f}")
print("(2) identical random laws on {-L..L}: min over laws of (d^2-Y^2)/J, worst excess Y^2-d^2 f(J/d^2)")
for L in (1,2,3):
  laws=[]
  while len(laws)<1500:
    w=rng.dirichlet(np.ones(L+1)*rng.choice([0.2,0.5,1,3]))
    p=np.concatenate([w[:0:-1]/2,[w[0]],w[1:]/2]); d=strip_law(p)
    if d>1e-3: laws.append((p,d))
  for J in (0.01,0.05,0.2,0.5,1,2):
    mr=1e9; ex=-1e9; n=0
    for p,d in laws:
        Y=Y_pair(p,p,J); 
        if Y>1e-6:
            mr=min(mr,(d*d-Y*Y)/J); n+=1
        ex=max(ex,Y*Y-d*d*f(J/(d*d)))
    print(f" L={L} J={J:5.2f}: min(d2-Y2)/J={mr:8.4f} (over {n} laws with Y>0)  max[Y2-d2 f]={ex:9.2e}")
print("(3) star: spin in C_d, N Ising, complex fields Re t>0: max over tests of ymax-sqrt(max(d^2-sum a^2,0))")
def star_poly(q,a,t):             # q ascending coefficients in w; returns polynomial in w
    N=len(a); out=np.zeros(len(q),dtype=complex)
    for e in np.ndindex(*(2,)*N):
        eps=2*np.array(e)-1; s=float(eps @ a); wt=np.exp(np.sum(t*eps))
        qq=q.copy(); # shift w->w+s
        coef=np.zeros(len(q),dtype=complex); pw=np.array([1.0+0j])
        for c in q:
            coef[:len(pw)]+=c*pw; pw=P.polymul(pw,[s,1.0])
        out+=wt*coef
    return out
worst=-1e9; tot=0
for it in range(6000):
    n=rng.integers(1,5); ds=rng.uniform(0,1,n); ds[0]=1; ts=rng.normal(0,rng.choice([.2,1,3]),n)
    q=np.array([1.0+0j])
    for t,d in zip(ts,ds): q=P.polymul(q,[(1j*t)**2-d*d,2j*t,1.0])
    N=rng.integers(1,4); a=rng.uniform(0,1,N)*rng.choice([.2,.5,1])
    t=rng.exponential(rng.choice([.03,.3,1.5]),N)+1j*rng.uniform(-3,3,N)
    r=P.polyroots(star_poly(q,a,t)); y=r.real.max(); Yb=math.sqrt(max(1-(a*a).sum(),0))
    worst=max(worst,y-Yb); tot+=1
print(f" tests={tot} worst excess={worst:.3e}")
print("(4) sharpness: partner Ising(1), first spin psi^m (spin one, strip a), a=J: limit arccosh(cosh a/cosh J); sqrt(a^2-J^2)")
for (q,J) in ((0.2,0.3),(0.2,0.7),(0.249,0.05)):
    a=math.acosh((1-2*q)/(2*q)); lim=math.acosh(math.cosh(a)/math.cosh(J))
    row=[]
    for m in (1,2,5,20,80):
        p=np.array([1.0]); 
        for _ in range(m): p=np.convolve(p,[q,1-2*q,q])
        p/=p.max()
        row.append(Y_pair(p,np.array([.5,0,.5]),J))   # partner on {-1,0,1} with mass only on +-1
    print(f" q={q} J={J} a={a:.4f} limit={lim:.4f} sqrt(a2-J2)={math.sqrt(a*a-J*J):.4f} Y(m=1,2,5,20,80)="+" ".join(f"{x:.4f}" for x in row))
