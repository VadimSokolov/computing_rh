"""Track E: rigorous moment tests without input zeros.
Run: timeout 600 python3 -I moment_tests.py
Third party dependency: python-flint. The /tmp/d1lib fallback contains only
an existing python-flint package; no external research scripts are imported.
"""
import sys, time, math, json
from pathlib import Path
sys.path.insert(0, '/tmp/e_lib')
sys.path.append('/tmp/d1lib')
from flint import arb, acb, arb_series, arb_mat, ctx
ctx.prec = 768
ctx.cap = 66
PI = arb.pi()
A = arb(4)
KAP = ((10-2*arb(5).sqrt()).sqrt()-2)/(arb(5).sqrt()-1)
START = time.monotonic()
LOG=[]
def say(*args):
    line=' '.join(str(x) for x in args)
    LOG.append(line)
    print(line, flush=True)
def out(x): return x.str(26)
def coeff(n): return [arb(0),arb(1),KAP,-KAP,arb(-1)][n%5]
def kernel(x, law):
    y=(2*x).exp()
    if law=='polya':
        return sum((4*(PI*n*n)**2*(arb(9)*x/2).exp()-6*PI*n*n*(arb(5)*x/2).exp())*(-PI*n*n*y).exp() for n in range(1,17))
    return 2*(arb(3)*x/2).exp()*sum(n*coeff(n)*(-PI*n*n*y/5).exp() for n in range(1,41))
def geo(c,N,power):
    n=arb(N+1)
    ratio=((n+1)/n)**power*(-c*(2*n+1)).exp()
    assert ratio < 1
    return n**power*(-c*n*n).exp()/(1-ratio)
def error(k,law,h=arb(0),b=arb(0)):
    # Bound both omitted theta terms on [0,A] and the integral on [A,infinity).
    # cosh(hx), |sinh(hx)| <= exp(hx), all h,b used are nonnegative.
    if law=='polya':
        inside=8*PI**2*(arb(9)*A/2+h*A+b*A*A).exp()*A**(k+1)/(k+1)*geo(PI,16,4)
        c=PI; exponent=arb(9)/2; C=16*PI**2
    else:
        inside=4*(arb(3)*A/2+h*A+b*A*A).exp()*A**(k+1)/(k+1)*geo(PI/5,40,1)
        c=PI/5; exponent=arb(3)/2; C=arb(8)
    # Log derivative is k/x+exponent+h+2*b*x-2*c*exp(2*x).
    # Its derivative is <=2*b-4*c*exp(2*A)<0 on x>=A.
    assert 2*b-4*c*(2*A).exp() < 0
    rate=2*c*(2*A).exp()-exponent-h-2*b*A-arb(k)/A
    assert rate > 0
    outside=C*A**k*(exponent*A+h*A+b*A*A-c*(2*A).exp()).exp()/rate
    return inside+outside

def integral_moment(k,law,h=arb(0),b=arb(0)):
    def fun(x,analytic):
        w=(h*x).cosh() if k%2==0 else (h*x).sinh()
        return 2*x**k*kernel(x,law)*w*(b*x*x).exp()
    val=acb.integral(fun,0,A,abs_tol=arb(2)**-680,rel_tol=arb(2)**-680,eval_limit=200000,depth_limit=25)
    assert val.is_finite() and val.imag.contains(0)
    return val.real+arb(0,error(k,law,h,b).upper())

def mgf_series(law,h=arb(0)):
    s=arb_series([arb(1)/2+h,1])
    if law=='polya':
        f=s*(s-1)/2*(-s*PI.log()/2).exp()*(s/2).gamma()*s.zeta()
    else:
        d=sum(coeff(r)*s.zeta(arb(r)/5) for r in range(1,5))*(-s*arb(5).log()).exp()
        f=(((s+1)/2)*(arb(5)/PI).log()).exp()*((s+1)/2).gamma()*d
    return f/f[0]
def p_from_moments(m):
    f=arb_series([m[k//2]/math.factorial(k) if k%2==0 else arb(0) for k in range(65)])
    g=f.log()
    return [arb(0)]+[(-1)**(j+1)*j*g[2*j] for j in range(1,33)]
def determinant(p,shift=11,size=4,divisor='uniform'):
    def q(n):
        den=arb(2*n).zeta() if divisor in ('uniform','sign') else arb(1)
        if divisor=='sign': den *= arb(2)**(2*n)-1
        return p[n]*arb(800 if divisor=='sign' else 200)**n/den
    mat=arb_mat([[q(shift+i+j) for j in range(size)] for i in range(size)])
    return mat.det(),mat
say('precision_bits',ctx.prec,'degree',64,'cutoff_x',A,'theta_terms_polya',16,'theta_terms_DH',40)
say('DH_kappa',out(KAP))
for law in ['polya','dh']:
    exact=mgf_series(law)
    assert all(exact[k].contains(0) for k in range(1,65,2))
    moments=[integral_moment(k,law) for k in range(0,65,2)]
    norm=moments[0]
    moments=[x/norm for x in moments]
    assert all(moments[j].overlaps(exact[2*j]*math.factorial(2*j)) for j in range(33))
    p=p_from_moments(moments)
    say(law,'normalization',out(norm),'variance',out(moments[1]),'kappa4',out(-12*p[2]),'kappa6',out(240*p[3]))
    for label in ['uniform','LY']:
        d,mat=determinant(p,divisor=label)
        say(law,label,'scaled_det_shift11_size4',out(d),'sign',('negative' if d<0 else 'positive' if d>0 else 'unresolved'))
        if law=='polya' and label=='uniform':
            assert d < 0
            say('polya_uniform_raw_det',out(d/arb(200)**56))
            for i in range(4): say('matrix_row',i,*[out(mat[i,j]) for j in range(4)])
    d,_=determinant(p,shift=10,size=12,divisor='sign')
    say(law,'sign_scaled_det_shift10_size12',out(d),'sign',('negative' if d<0 else 'positive' if d>0 else 'unresolved'))
    if law=='polya': assert d<0
    d,_=determinant(p,shift=3,size=6,divisor='uniform')
    say(law,'uniform_scaled_det_shift3_size6',out(d),'sign',('negative' if d<0 else 'positive' if d>0 else 'unresolved'))
    if law=='polya':
        assert all(p[j]>0 for j in range(1,33))
        say(law,'alternating_cumulants_through_order',64,'all_pass')
    # GHS is stronger than the one variable Lee and Yang property.
    # The tilted third cumulant is a necessary GS inequality.
    for b,h in [(arb(0),arb(1)),(arb(0),arb(5)),(arb(1),arb(1))]:
        v=[integral_moment(k,law,h,b) for k in range(4)]
        c3=v[3]/v[0]-3*v[2]*v[1]/v[0]**2+2*(v[1]/v[0])**3
        say(law,'GHS_b_h',b,h,'third_cumulant',out(c3),'pass' if c3<0 else 'FAIL' if c3>0 else 'unresolved')
    for j in [1,2,3,11,12,13,14,15,16,17,32]: say(law,'moment',2*j,out(moments[j]))
say('elapsed_seconds',round(time.monotonic()-START,3))
Path('results').mkdir(exist_ok=True)
Path('results/numbers.txt').write_text('\n'.join(LOG)+'\n')
