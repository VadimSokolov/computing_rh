import sys
sys.path.insert(0, '/tmp/e_lib')
sys.path.append('/tmp/d1lib')
import flint
from flint import arb, arb_series, arb_mat, acb, ctx
print(sys.version)
print('flint', flint.__version__, flint.__file__)
print(arb_series.zeta.__doc__)
print(acb.integral.__doc__[:2400])
ctx.prec=512
ctx.cap=66
s=arb_series([arb(1)/2,1])
f=s*(s-1)/2*(-s*arb.pi().log()/2).exp()*(s/2).gamma()*s.zeta()
f=f/f[0]
g=f.log()
p=[arb(0)]+[(-1)**(j+1)*j*g[2*j] for j in range(1,33)]
for label, q in [('LY',p),('uniform',[arb(0)]+[p[j]/arb(2*j).zeta() for j in range(1,33)]),('sign',[arb(0)]+[p[j]/((arb(2)**(2*j)-1)*arb(2*j).zeta()) for j in range(1,33)])]:
    found=False
    for size in range(2,9):
        for shift in range(1,33-2*(size-1)):
            mat=arb_mat([[q[shift+i+j] for j in range(size)] for i in range(size)])
            det=mat.det()
            if det < 0:
                print('VIOLATION',label,'size',size,'shift',shift,det)
                found=True
                break
        if found: break
    print(label,'done', 'q1-3',q[1:4])
print('kappa4',24*g[4])
