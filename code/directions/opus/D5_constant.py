# D5 experiment 5: how sharp is the lattice constant M*kappa?  Partner on {-1,0,1} with P(+-1)=q0 (Lee and Yang for q0>=1/4), field h2=0 (limit y2->0+).
# First spin local: (w^2-1)^m.  F_m(w)=(1-2q0)Q(w)+q0 Q(w+J)+q0 Q(w-J).  Theorem 2 gives y^2<=1-J^2/2; variance would give 1-2 q0 J^2; two point limit 1-J^2.
import mpmath as mp
mp.mp.dps=120
def polymul(a,b):
    c=[mp.mpf(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]+=x*y
    return c
def shift(p,s):                 # p(w+s), ascending
    out=[mp.mpf(0)]*len(p); pw=[mp.mpf(1)]
    for c in p:
        for i,x in enumerate(pw): out[i]+=c*x
        pw=polymul(pw,[mp.mpf(s),mp.mpf(1)])
    return out
def ymax(m,q0,J):
    Q=[mp.mpf(1)]
    for _ in range(m): Q=polymul(Q,[mp.mpf(-1),mp.mpf(0),mp.mpf(1)])
    A=shift(Q,J);B=shift(Q,-J)
    F=[(1-2*q0)*Q[i]+q0*A[i]+q0*B[i] for i in range(len(Q))]
    r=mp.polyroots(F[::-1],maxsteps=5000,extraprec=800)
    return max(mp.re(x) for x in r)
print("q0     J    1-J^2/2  1-2q0J^2  1-J^2 | y^2 for m=1,2,3,4,6,8,12,16,24")
for q0 in (mp.mpf(1)/4,mp.mpf(3)/10,mp.mpf(3)/8,mp.mpf(45)/100,mp.mpf(49)/100):
    for J in (mp.mpf(1)/5,mp.mpf(1)/2):
        row=[ymax(m,q0,J) for m in (1,2,3,4,6,8,12,16,24)]
        print(f"{float(q0):5.3f} {float(J):4.2f} {float(1-J*J/2):8.4f} {float(1-2*q0*J*J):8.4f} {float(1-J*J):8.4f} | "+" ".join(f"{float(max(y,0))**2:.4f}" for y in row))
