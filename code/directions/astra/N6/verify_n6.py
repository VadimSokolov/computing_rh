"""Independent checks of cosh_spin.tex. No input is read from the note folder.

Run with timeout 600 python3 -I verify_n6.py MODE.
Ball enclosures are used for moments and signs. Sampled winding numbers and
zero scans are explicitly numerical evidence, not completeness certificates.
"""
import sys
sys.path.insert(0, '/tmp/cs_lib')
import json
import math
import time
from pathlib import Path
import mpmath as mp
from flint import arb, acb, arb_poly, ctx

ctx.prec = 192
mp.mp.dps = 65
PI = arb.pi()
OUT = Path(__file__).resolve().parent / 'results'
OUT.mkdir(exist_ok=True)

def record(name, data):
    data['seconds'] = time.monotonic() - START
    (OUT / (name + '.json')).write_text(json.dumps(data, indent=2) + '\n')
    print(json.dumps(data, indent=2), flush=True)

def val(x):
    return float(x.mid())

def integrate(fun, lo, hi):
    ans = acb.integral(lambda z, analytic: fun(z), lo, hi,
                       rel_tol=arb(2)**-150, abs_tol=arb(2)**-165)
    assert ans.is_finite(), ans
    return ans.real

def cosh_weight(x, a=1, b=0):
    return (-a * (x.cosh() - 1) - b*x*x).exp()

def star_weight(u):
    return (2*PI**2*(9*u).cosh()-3*PI*(5*u).cosh()) * (-2*PI*(4*u).cosh()).exp()

def star_tail(k, h=0, L=2):
    slope = 4*PI*(arb(4)*L).exp() - (k+9+h)
    assert slope > 0
    return 2*PI**2*((k+9+h)*L-PI*(arb(4)*L).exp()).exp()/slope

def moments_star(kmax=10, h=0):
    raw=[]
    for k in range(kmax+1):
        if h == 0 and k % 2:
            raw.append(arb(0))
            continue
        fn=lambda u: u**k * star_weight(u) * ((h*u).exp()+(-1)**k*(-h*u).exp())
        r=integrate(fn, 0, 2)
        tail=2*star_tail(k, abs(h))
        r += arb(0, tail)
        raw.append(r)
    return [r/raw[0] for r in raw]

def cumulants(m):
    k=[arb(0)]
    for n in range(1,len(m)):
        k.append(m[n]-sum(math.comb(n-1,j-1)*k[j]*m[n-j] for j in range(1,n)))
    return k

def F(t, n=1):
    t=acb(t)
    z=acb(2*PI*n*n)
    q=acb(0,1)*t/2
    return 2*PI**2*n**4*z.bessel_k(acb(9)/4+q)-3*PI*n*n*z.bessel_k(acb(5)/4+q)

def transform(t, N=1):
    t=acb(t)
    return sum((F(t,n)+F(-t,n))/4 for n in range(1,N+1))

def real_transform(t, N=1):
    return sum(F(t,n).real/2 for n in range(1,N+1))

def theta(t):
    return mp.siegeltheta(t)

def gamma_infty(t):
    s=acb(arb(1)/2,arb(t))
    return s*(s-1)/2 * (-s/2*PI.log()).exp() * (s/2).gamma()

def sign(f):
    if f>0: return 1
    if f<0: return -1
    raise ValueError('Indeterminate sign: '+str(f))

def bisect(fun, a,b, iterations=44):
    a=arb(a); b=arb(b); sa=sign(fun(a))
    assert sa!=sign(fun(b))
    for _ in range(iterations):
        c=(a+b)/2
        if sign(fun(c))==sa: a=c
        else: b=c
    return [val(a),val(b)]

def scan(fun, a,b, step):
    brackets=[]
    last=a; slast=sign(fun(a))
    for j in range(1,round((b-a)/step)+1):
        t=a+j*step; st=sign(fun(t))
        if st!=slast: brackets.append([last,t])
        slast=st; last=t
    return brackets

def folded(t, N):
    t=acb(t)
    ans=acb(0)
    for n in range(1,N+1):
        a=acb(PI*n*n)
        for c,b in [(2*PI**2*n**4,arb(9)/4),(-3*PI*n*n,arb(5)/4)]:
            for sgn in [-1,1]:
                nu=b+sgn*acb(0,1)*t/2
                ans += c/4 * (-nu*a.log()).exp()*a.gamma_upper(nu)
    return ans.real

def dh_kappa():
    return ((10-2*arb(5).sqrt()).sqrt()-2)/(arb(5).sqrt()-1)

def dh(s):
    s=acb(s); k=dh_kappa()
    coeff=[arb(1),k,-k,arb(-1)]
    return (-s*arb(5).log()).exp()*sum(coeff[r-1]*s.zeta(acb(r)/5) for r in range(1,5))

def do_moments():
    data={}
    m=moments_star(); k=cumulants(m)
    data['star_cumulants']={str(j):str(k[j]) for j in range(2,11,2)}
    data['star_excess']=str(k[4]/k[2]**2)
    data['star_moments']={str(j):str(m[j]) for j in range(0,11,2)}
    raw=[]
    for j in range(7):
        if j%2: raw.append(arb(0)); continue
        r=2*integrate(lambda x:x**j*cosh_weight(x,1,10),0,8)
        tail=2*(arb(1)+j*8-arb(8).exp()/2).exp()/(arb(8).exp()/2-j)
        raw.append(r+arb(0,tail))
    kd=cumulants([r/raw[0] for r in raw])
    assert kd[6]<0
    data['damped_cosh']={str(j):str(kd[j]) for j in [2,4,6]}
    # Algebraic numerator for the global GHS inequality, derived afresh.
    x=arb_poly([0,1]); r=3/(2*PI)
    q=x**9-r*x**7-r*x**2+1
    e=lambda p:x*p.derivative()
    q1=e(q); q2=e(q1); q3=e(q2)
    h=q3*q*q-3*q2*q1*q+2*q1**3
    p=8*PI*(x**4-1)*q**3-x*x*h
    shifted=[sum(p[i]*math.comb(i,j) for i in range(j,len(p))) for j in range(1,32)]
    assert all(c>0 for c in shifted)
    data['ghs_shifted_coefficient_lower_integers']=[math.floor(float(c.lower())) for c in shifted]
    assert all(c>b for c,b in zip(shifted,data['ghs_shifted_coefficient_lower_integers']))
    data['ghs_coefficient_balls']=[str(c) for c in shifted]
    upp=lambda u:8*p((2*arb(u)).exp())/((4*arb(u)).exp()*q((2*arb(u)).exp())**3)
    data['U_third']={u:str(upp(u)) for u in ['0.00000001','0.0001','0.01','0.02','2']}
    data['bessel_quadrature_checks']={}
    for xarg in [0,10]:
        direct=2*integrate(lambda u:star_weight(u)*(xarg*u).cos(),0,2)+arb(0,2*star_tail(0))
        exact=real_transform(arb(xarg)/2)
        data['bessel_quadrature_checks'][str(xarg)]={'quadrature':str(direct),'bessel':str(exact),'difference':str(direct-exact)}
    data['Kthird_grid']={}
    vals=[]
    for j in range(1,241):
        hval=arb(j)/2
        cur=cumulants(moments_star(3,hval))[3]
        assert cur<0
        vals.append((val(cur),float(hval)))
    data['Kthird_grid']['maximum']=max(vals)
    data['Kthird_grid']['at_200']=str(cumulants(moments_star(3,200))[3])
    data['harmonic_cumulants']={str(j):str(arb([-2,16,-272,7936][j//2-2])*arb(j).zeta()) for j in [4,6,8,10]}
    record('moments',data)

def do_zeros():
    data={}
    for N in [1,2,3]:
        limit=700 if N<3 else 100
        br=scan(lambda t:real_transform(t,N),0,limit,.25)
        data['N'+str(N)]={'sign_changes_to_100':sum(b[1]<=100 for b in br),
                         'sign_changes_to_260':sum(b[1]<=260 for b in br),
                         'sign_changes_200_700':sum(b[0]>=200 for b in br),
                         'brackets':br}
        small=[bisect(lambda t:real_transform(t,N),*b) for b in br[:8]]
        data['N'+str(N)]['first_eight_brackets']=small
        if N==1:
            idx=[0,5]+[i for i,b in enumerate(br) if 257<b[0]<260]
            phases=[]
            for i in idx:
                box=bisect(lambda t:real_transform(t,1),*br[i]); t=sum(box)/2
                p=theta(mp.mpf(t))/mp.pi
                phases.append({'t':t,'theta_over_pi':str(p),'fraction':str(mp.frac(p)),'t_deficit':str(t*(mp.mpf('.5')-mp.frac(p)))})
            data['phase_checks']=phases
    main=lambda t:arb(t)/(2*PI)*(arb(t)/(2*PI)).log()-arb(t)/(2*PI)+arb(7)/8
    data['main_counts']={str(t):str(main(t)) for t in [100,260]}
    data['main_counts']['200_to_700']=str(main(700)-main(200))
    data['error_E']={}
    for t in [17.846,18.5,30,36.64,100,700]:
        err=2*F(t)/gamma_infty(t)-1
        data['error_E'][str(t)]={'abs':str(abs(err)),'t_abs':str(t*abs(err)),'phase':str((1+err).arg())}
    data['gram_intervals_completed_to700']=int(mp.floor(theta(700)/mp.pi))
    data['theta700_overpi']=str(theta(700)/mp.pi)
    data['ratios_of_maxima']={}
    for a,b in [(20,40),(60,80),(100,120),(200,220),(680,700)]:
        grid=[arb(a)+arb(j)/10 for j in range(10*(b-a)+1)]
        max1=max(val(abs(F(t,1).real)) for t in grid)
        max2=max(val(abs(F(t,2).real)) for t in grid)
        data['ratios_of_maxima'][str((a,b))]=max2/max1
    record('zeros',data)

def do_winding():
    data={}
    for N in [1,2]:
        for step in [.5,.25]:
            pts=[]
            corners=[complex(200,-1),complex(700,-1),complex(700,1),complex(200,1),complex(200,-1)]
            for a,b in zip(corners,corners[1:]):
                n=round(abs(b-a)/step)
                pts.extend(a+(b-a)*j/n for j in range(n))
            pts.append(pts[0])
            phases=[]
            uncertain=0
            for z in pts:
                zz=transform(acb(z.real,z.imag),N)
                if abs(zz).contains(0): uncertain+=1
                # The midpoint of an argument ball crossing its branch cut is
                # not the argument of the midpoint. Use the value's midpoint.
                scale=abs(zz)
                phases.append(math.atan2(val(zz.imag/scale),val(zz.real/scale)))
            increments=[math.remainder(b-a,2*math.pi) for a,b in zip(phases,phases[1:])]
            data[str((N,step))]={'sampled_winding':math.fsum(increments)/(2*math.pi),
                                'max_sample_phase_increment':max(abs(d) for d in increments),'evaluations':len(pts),'uncertain_values':uncertain}
    record('winding',data)

def do_folded():
    data={}
    for N in [1,2,3,4]:
        br=scan(lambda t:folded(t,N),0,200,.25)
        cusp=sum(PI*n*n*(-PI*n*n).exp()*(30*PI*n*n-15-8*PI**2*n**4) for n in range(1,N+1))
        data[str(N)]={'positive_real_sign_changes_to200':len(br),
                      'last_bracket':bisect(lambda t:folded(t,N),*br[-1]),
                      'cusp_slope':str(cusp),'value_200':str(folded(200,N)),'brackets':br}
    record('folded',data)

def do_dh():
    k=mp.mpf(str(dh_kappa().mid()).split(' ')[0].strip('[]'))
    def d(s):
        return mp.power(5,-s)*sum(c*mp.zeta(s,mp.mpf(r)/5) for r,c in [(1,1),(2,k),(3,-k),(4,-1)])
    root=mp.findroot(d,mp.mpc('.8085','85.6993'))
    sb=acb(str(root.real),str(root.imag))
    data={'kappa':str(dh_kappa()),'root_computed':str(root),'ball_value_at_rounded_root':str(dh(sb)),
          'real_axis_positive_kernel_bound':'tail / first <= exp(-3*pi/5)*(2+5/(2*pi)) < 1'}
    data['tail_bound_at1']=str((-3*PI/5).exp()*(2+5/(2*PI)))
    # Moments for U with Fourier variable t: H(u)=exp(3u/2)*theta_D(exp(2u)).
    # On u>=0, truncate n<=12. The omitted terms are bounded analytically
    # by sum n exp(-pi*n*n/5) times exp(-positive u); this check is numerical.
    coeff=lambda n:[arb(0),arb(1),dh_kappa(),-dh_kappa(),arb(-1)][n%5]
    def density(u):
        return (3*u/2).exp()*sum(coeff(n)*n*(-PI*n*n/5*(2*u).exp()).exp() for n in range(1,13))
    c=PI/5
    missing_n=2/(13*c)*(-169*c).exp()*(1+1/(26*c))
    raws=[]
    for j in range(11):
        if j%2:
            raws.append(arb(0)); continue
        tail_u=2*arb('1.425')*((j+arb(3)/2)*4-c*arb(8).exp()).exp()/(2*c*arb(8).exp()-(j+arb(3)/2))
        raws.append(2*integrate(lambda u:u**j*density(u),0,4)+arb(0,missing_n+tail_u))
    kd=cumulants([v/raws[0] for v in raws])
    data['dh_cumulants_with_tail_balls']={str(j):str(kd[j]) for j in range(2,11,2)}
    data['dh_excess_with_tail_ball']=str(kd[4]/kd[2]**2)
    data['moment_tail_bound']=str(missing_n+tail_u)
    for j in range(2,11,2): assert (-1)**(j//2-1)*kd[j]>0
    record('dh',data)

def do_summary():
    for p in sorted(OUT.glob('*.json')):
        if p.stem=='boundary': continue
        d=json.loads(p.read_text())
        def trim(v):
            if isinstance(v,dict):
                return {k:trim(x) for k,x in v.items() if k not in ['brackets','ghs_coefficient_balls','star_moments','gram_interval_assignments_ball_checked']}
            return v
        print(p.name,json.dumps(trim(d),indent=2))

def do_extra():
    # Check the below-100 root-shift assertion at the last ten roots.
    d=json.loads((OUT/'zeros.json').read_text())
    roots={}
    for N in [1,2]:
        br=[b for b in d['N'+str(N)]['brackets'] if b[1]<=100]
        roots[N]=[sum(bisect(lambda t:real_transform(t,N),*b))/2 for b in br]
    out={'root_shifts_below100':[{'index':j+1,'N1':a,'N2':b,'shift':b-a} for j,(a,b) in enumerate(zip(roots[1],roots[2])) if j>=19]}
    def xi(s):
        return s*(s-1)/2*mp.power(mp.pi,-s/2)*mp.gamma(s/2)*mp.zeta(s)
    kk=[mp.diff(lambda h:mp.log(xi(mp.mpf('.5')+h/2)),0,j) for j in [2,4]]
    out['polya_cumulants_computed']={str(j):str(k) for j,k in zip([2,4],kk)}
    out['polya_excess_computed']=str(kk[1]/kk[0]**2)
    out['first_eight_zeta_ordinates']=[str(acb.zeta_zero(j).imag) for j in range(1,9)]
    # Gram endpoint signs, ball values at high precision approximations to g_k.
    # Their tiny input enclosure is not a validated enclosure of the Gram root.
    bad=[]
    for j in range(414):
        g=mp.grampoint(j)
        v=real_transform(arb(str(g)))
        if sign(v)!=(-1)**(j+1): bad.append(j)
    out['gram_endpoint_sign_mismatches_to_g413']=bad
    out['g413_computed']=str(mp.grampoint(413))
    out['g414_computed']=str(mp.grampoint(414))
    def theta_ball(t):
        t=arb(t)
        return acb(arb(1)/4,t/2).lgamma().imag-t/2*PI.log()
    assigned={}
    for lo,hi in d['N1']['brackets']:
        if lo<17: continue
        left=theta_ball(lo)/PI; right=theta_ball(hi)/PI
        k=math.floor(val(left))
        assert left>k and right<k+1
        assigned[k]=assigned.get(k,0)+1
    assert all(assigned[k]==1 for k in range(414))
    out['gram_interval_assignments_ball_checked']=assigned
    derivative=(acb(arb(1)/4,arb(17)/2).digamma().real-PI.log())/2
    assert derivative>0
    out['theta_prime_17']=str(derivative)
    record('extra',out)

def do_boundary():
    out={}
    for N in [1,2]:
        out[str(N)]=[]
        for x in [200,250,400,600,700]:
            for y in [-1,0,1]:
                z=acb(x,y); f=transform(z,N)
                out[str(N)].append({'z':[x,y],'f':str(f),'arg':str(f.arg()),'abs':str(abs(f)),
                                    'real_check':str(f.real-real_transform(x,N)) if y==0 else None})
    record('boundary',out)

def do_contour(low=False):
    out={}
    cases=[(1,1,700)] if low else [(1,200,700),(2,200,700)]
    for N,left,right in cases:
        evaluations=[0]
        def fn(z):
            evaluations[0]+=1
            if low and z.real < 20:
                ans=2*acb.integral(lambda u,analytic:star_weight(u)*(2*z*u).cos(),0,2,
                                   rel_tol=arb(2)**-130,abs_tol=arb(2)**-145)
                tail=2*star_tail(0,2*abs(z.imag).upper())
                return ans+acb(arb(0,tail),arb(0,tail))
            return transform(z,N)
        def piece(a,b,fa,fb,depth=0):
            mid=(a+b)/2
            box=acb(arb(mid.real,abs(b.real-a.real)/2),arb(mid.imag,abs(b.imag-a.imag)/2))
            cover=fn(box)
            phase=(fb/fa).arg()
            # A convex complex ball avoiding zero contracts the whole image
            # segment to its chord without changing the winding number.
            if cover.is_finite() and not cover.contains(0) and abs(phase)<PI:
                return phase,1,depth
            if depth>=14:
                raise RuntimeError('Contour subdivision limit '+str((a,b,cover)))
            fm=fn(acb(mid.real,mid.imag))
            p1,n1,d1=piece(a,mid,fa,fm,depth+1)
            p2,n2,d2=piece(mid,b,fm,fb,depth+1)
            return p1+p2,n1+n2,max(d1,d2)
        corners=[complex(left,-1),complex(right,-1),complex(right,1),complex(left,1),complex(left,-1)]
        total=arb(0); pieces=0; maxdepth=0
        for a,b in zip(corners,corners[1:]):
            count=round(abs(b-a)/.5)
            fa=fn(acb(a.real,a.imag))
            for j in range(count):
                lo=a+(b-a)*j/count; hi=a+(b-a)*(j+1)/count
                fb=fn(acb(hi.real,hi.imag))
                p,num,depth=piece(lo,hi,fa,fb)
                total+=p; pieces+=num; maxdepth=max(maxdepth,depth); fa=fb
        winding=total/(2*PI)
        nearest=round(val(winding))
        assert winding.contains(nearest) and abs(winding-nearest)<arb(1)/4
        out[str(N)]={'winding_ball':str(winding),'count':nearest,'pieces':pieces,
                     'max_subdivision_depth':maxdepth,'evaluations':evaluations[0],
                     'real_bounds':[left,right]}
        print('Completed certified contour',N,nearest,flush=True)
    record('contour_low' if low else 'contour',out)

START=time.monotonic()
if __name__=='__main__':
    mode=sys.argv[1]
    {'moments':do_moments,'zeros':do_zeros,'winding':do_winding,'folded':do_folded,'dh':do_dh,'summary':do_summary,'extra':do_extra,'boundary':do_boundary,'contour':do_contour,'contour_low':lambda:do_contour(True)}[mode]()
