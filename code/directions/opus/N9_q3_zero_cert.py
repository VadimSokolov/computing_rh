# Q3: the note's grid test as a corollary of the zeros: T_m >= T^(K)_m - tau I, certified in balls from Arb's zeros.
import sys; sys.path.append('./misc/research/directions/opus')
from N9_common import *
t0=time.time(); ctx.prec=128
K=int(sys.argv[1]) if len(sys.argv)>1 else 2500
zs=acb.zeta_zeros(1,K); gam=[z.imag for z in zs]; print('Arb zeros: K=%d, gamma_K=%.3f, %.1fs; on the line: %s'%(K,fl(gam[-1]),time.time()-t0,all(abs(fl(z.real)-0.5)<1e-30 for z in zs)))
d=arb(1)/40; df=0.025; M=482
gf=np.array([fl(g) for g in gam])
# build c^(K) in balls with the angle addition recurrence avoided: direct cos
c=[arb(0) for _ in range(M)]
for g in gam:
    w=4*(g*d/2).sin()**2/(g*g); th=g*d
    for m in range(M):
        c[m]+=w*(arb(m)*th).cos()
print('c^(K) built, max radius %.2e, %.1fs'%(max(fl(x.rad()) for x in c),time.time()-t0))
H=3e12
print(' a | m | lam_min(T^(K)) float | tau(H=3e12, beta<=1/2) | certified T^(K)-tau I > 0 | info')
cf=np.array([fl(x) for x in c])
for a in [6,12.05]:
    m=int(round(a/df)); lam=lam_min(cf,m)
    tau=math.sinh(a/2)/math.tanh(df/4)*S_tail(H)
    ok,info=chol_cert(c[:m],float(tau))
    print(' %g | %d | %.5e | %.5e | %s | %s'%(a,m,lam,tau,ok,['%.2e'%x for x in info]))
print('total %.1f'%(time.time()-t0))
