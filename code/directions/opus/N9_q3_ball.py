# Q3: the discretised test in balls, delta=1/40, against the certified ladder.
import sys; sys.path.append('./misc/research/directions/opus')
from N9_common import *
t0=time.time()
ctx.prec=256; C=consts(); pps=pp_ball(250000)
delta=arb(1)/40; M=482; d=0.025
P=Psi_grid_ball(delta,M,pps,C); c=cov_from_Psi_ball(P)
print('c at 256 bits, max radius %.2e'%max(fl(x.rad()) for x in c))
cf=np.array([fl(x) for x in c])
G=load_zeros(8000); print('N(pi/delta)=N(%.2f)=%d, gamma_241=%.2f'%(math.pi/d,int(np.sum(G<math.pi/d)),G[240]))
D1={2:(1.31651e-3,1.329433e-3,7.19318e-2,7.372781e-2),3:(5.51026e-8,5.535508e-8,1.48233e-5,1.528437e-5),4:(7.27636e-13,7.559845e-13,4.52355e-10,5.012206e-10),5:(8.71617e-18,9.344959e-18,9.03483e-15,1.030676e-14),7:(5.72322e-28,7.135123e-28,1.35675e-24,1.610952e-24),8:(2.67372e-33,4.373503e-33,9.85504e-30,1.245146e-29),9:(1.97280e-38,2.765132e-38,7.98786e-35,1.026469e-34)}
rows=[('log %d'%X,math.log(X),X) for X in D1]+[('%g'%a,a,None) for a in [3,4,5,6,8,10,12,12.05]]
print('\n a | X=e^a | m | certified lower < lam_min(T_m) < ball Rayleigh upper | 2 lam_min/delta | parity, |u_end|, peak freq | ladder margin (D1) | digits 4piX/ln10')
for lab,a,X in rows:
    m=int(math.floor(a/d+1e-9))
    lam,u=lam_min(cf,m,want_vec=True)
    ok,info=chol_cert(c[:m],lam*(1-1e-4))
    up=rayleigh_ball(c[:m],u)
    par='even' if u[0]*u[-1]>0 else 'odd'
    Xv=math.exp(a)
    lad=('%.3g (even), %.3g (odd)'%(D1[X][0],D1[X][2])) if X in D1 else 'none certified'
    print('%s | %.4g | %d | %.6e %s < %.6e | %.4e | %s %.2e %.2f | %s | %.4g'%(lab,Xv,m,lam*(1-1e-4),'CERT' if ok else 'FAIL',ub(up),2*lam/d,par,abs(u[0]),peak_freq(u,d),lad,4*math.pi*Xv/math.log(10)))
print('time %.1f'%(time.time()-t0))
# Levinson in balls: needs high precision because of wrapping
ctx.prec=6000; C=consts(); pps=pp_ball(250000)
delta=arb(1)/40
P=Psi_grid_ball(delta,M,pps,C); ch=cov_from_Psi_ball(P)
ks,vs,fail=levinson_ball(ch)
print('\nLevinson ball 6000 bits: first failure:',fail,' radius of last v %.2e, last kappa %.2e; time %.1f'%(fl(vs[-1].rad()),fl(ks[-1].rad()),time.time()-t0))
akap=[abs(fl(k)) for k in ks]; i=int(np.argmax(akap))
print('max|kappa| = %.8f at k=%d; c(0)=Psi(0.025)=%.8f'%(akap[i],i,fl(ch[0])))
order=np.argsort(akap)[::-1][:14]; print('top |kappa|:',[(int(k),round(float(k*d),3),round(fl(ks[k]),4)) for k in order])
vv=[fl(v) for v in vs]; print('min v/v0 = %.8f at k=%d; all v>0 certified: %s; v_last=%.6e'%(min(vv)/vv[0],int(np.argmin(vv)),all(v>0 for v in vs),vv[-1]))
print('max|kappa_k| over k>=20: %.4f; over k>=200: %.4f, last k=%d'%(max(akap[20:]),max(akap[200:] or [0]),len(akap)-1))
e241=E_law(241,6.0); print('\nD2 law at the same dimension: N=241 even cos modes, L=6 (a=12): E=%.1f nats = 10^-%.1f, b=%.1f, phi=%.3f'%(e241[0],e241[0]/math.log(10),e241[1],e241[2]))
X12=math.exp(12); print('rung a=12: X=%.1f, 4 pi X=%.1f nats=10^-%.0f, N*=2X log X=%.3e modes, digits 5.46X=%.3e, bits 18.13X=%.3e'%(X12,4*math.pi*X12,4*math.pi*X12/math.log(10),2*X12*12,5.46*X12,18.13*X12))
print('total %.1fs'%(time.time()-t0))
