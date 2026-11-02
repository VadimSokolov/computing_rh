import sys
sys.path.insert(0,'/tmp/e_lib')
sys.path.append('/tmp/d1lib')
from flint import arb, arb_series, arb_mat,ctx
ctx.prec=2048
ctx.cap=130
s=arb_series([arb(1)/2,1])
f=s*(s-1)/2*(-s*arb.pi().log()/2).exp()*(s/2).gamma()*s.zeta()
g=(f/f[0]).log()
p=[arb(0)]+[(-1)**(j+1)*j*g[2*j] for j in range(1,65)]
for label in ['uniform','sign']:
    q=[arb(0)]+[p[j]*arb(800 if label=='sign' else 200)**j/(arb(2*j).zeta()*((arb(2)**(2*j)-1) if label=='sign' else 1)) for j in range(1,65)]
    found=False
    for largest in range(3,65):
        for size in range(2,min(12,(largest+1)//2)+1):
            shift=largest-2*(size-1)
            d=arb_mat([[q[shift+i+j] for j in range(size)] for i in range(size)]).det()
            if d<0:
                print(label,'largest',largest,'size',size,'shift',shift,'det',d.str(30),flush=True)
                found=True
                break
        if found:break
    if not found:print(label,'none resolved',flush=True)
