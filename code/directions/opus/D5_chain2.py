# D5 experiment 9: (1) odd chains in the local model, large J limit against u_m=(m+1)/(2^(m+3)-3m-7), N=2m+1.
# (2) pairs (a,b) in the local model, large J limit against the largest zero of d^(2b)(w^2-1)^a.
# (3) law level closure coupling of the pairs (a,b) of powers of spin one. mpmath, not certified.
import mpmath as mp
mp.mp.dps=80
def pmul(a,b):
    c=[mp.mpf(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]+=x*y
    return c
def padd(a,b):
    n=max(len(a),len(b)); c=[mp.mpf(0)]*n
    for i,x in enumerate(a): c[i]+=x
    for i,x in enumerate(b): c[i]+=x
    return c
def der(p,k=1):
    for _ in range(k):
        p=[i*p[i] for i in range(1,len(p))] or [mp.mpf(0)]
    return p
def rootsY(c):
    while len(c)>1 and abs(c[-1])<mp.mpf(10)**(-70)*max(abs(x) for x in c): c=c[:-1]
    r=mp.polyroots(c[::-1],maxsteps=8000,extraprec=1500)
    return r
Q=[mp.mpf(-1),mp.mpf(0),mp.mpf(1)]
Qd=[Q,der(Q),der(Q,2)]
def chain(N,J):
    A={k:[x*J**k/mp.factorial(k) for x in Qd[k]] for k in range(3)}
    for i in range(2,N):
        B={}
        for k in range(3):
            s=[mp.mpf(0)]
            for kp in range(3):
                if kp+k<=2: s=padd(s,pmul(A[kp],Qd[kp+k]))
            B[k]=[x*J**k/mp.factorial(k) for x in s]
        A=B
    s=[mp.mpf(0)]
    for kp in range(3): s=padd(s,pmul(A[kp],Qd[kp]))
    return s
print("(1) odd chains, largest real part squared of a real zero or any zero, local model delta=1")
for m in (1,2,3,4):
    N=2*m+1; lim=mp.mpf(m+1)/(2**(m+3)-3*m-7)
    vals=[]
    for J in (mp.mpf(1),mp.mpf(10),mp.mpf(100),mp.mpf(10000)):
        r=rootsY(chain(N,J)); vals.append(max(mp.re(z)**2 for z in r))
    print(f"N={N} u_m={mp.nstr(lim,8)} y^2 at J=1,10,100,1e4:",[mp.nstr(v,7) for v in vals])
print("(2) pairs (a,b), local model")
def pow1(a):
    p=[mp.mpf(1)]
    for _ in range(a): p=pmul(p,Q)
    return p
def Zpair(a,b,J):
    P1,P2=pow1(a),pow1(b); s=[mp.mpf(0)]; k=0
    d1,d2=P1,P2; t=mp.mpf(1)
    while True:
        s=padd(s,[x*t for x in pmul(d1,d2)])
        k+=1
        if k>min(2*a,2*b): break
        d1,d2=der(d1),der(d2); t=t*J/k
    return s
for (a,b) in ((2,1),(3,1),(4,1),(3,2),(4,2),(5,2),(4,3),(5,3)):
    top=der(pow1(a),2*b); rr=mp.polyroots(top[::-1],maxsteps=4000,extraprec=800) if len(top)>1 else []
    lim=max(mp.re(z)**2 for z in rr) if rr else mp.mpf(0)
    vals=[max(mp.re(z)**2 for z in rootsY(Zpair(a,b,mp.mpf(J)))) for J in (0.1,1,10,1000)]
    print(f"({a},{b}) limit {mp.nstr(lim,7)} y^2 at J=.1,1,10,1000:",[mp.nstr(v,7) for v in vals])
mp.mp.dps=50
def Ylaw(a,b,s,Jp):
    q=1/(2*(1+mp.cosh(s))); p=[q,1-2*q,q]
    def pw(m):
        out=[mp.mpf(1)]
        for _ in range(m): out=pmul(out,p)
        return out
    P1,P2=pw(a),pw(b); c=[mp.mpf(0)]*(len(P1)+len(P2)-1)
    for i,x in enumerate(P1):
        for j,y in enumerate(P2): c[i+j]+=x*y*mp.exp(Jp*(i-a)*(j-b))
    M=max(abs(x) for x in c); c=[x/M for x in c]
    r=mp.polyroots(c[::-1],maxsteps=4000,extraprec=600)
    return max(abs(mp.log(abs(z))) for z in r)
print("(3) law closure J_phys (lattice spacing 1, strip s); ratio to s and to s^2; bisection, monotonicity assumed")
for (a,b) in ((1,1),(2,2),(2,1),(3,1)):
    for s in (mp.mpf(1),mp.mpf('.3'),mp.mpf('.1'),mp.mpf('.03')):
        lo,hi=mp.mpf(0),mp.mpf(8)
        for _ in range(26):
            mid=(lo+hi)/2
            if Ylaw(a,b,s,mid)>mp.mpf(10)**(-14): lo=mid
            else: hi=mid
        print(f"({a},{b}) s={float(s):5.2f} J*={float(hi):9.5f} J*/s={float(hi/s):8.4f} J*/s^2={float(hi/s/s):9.3f}")
