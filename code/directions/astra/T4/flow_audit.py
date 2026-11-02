"""Independent finite-flow and Davenport-Heilbronn audit.

Run with PYTHONPATH=./pylib python3 flow_audit.py MODE, MODE in
ball, table, dh, dhflow, zeros, all. Bounded to 600 seconds per invocation.
Only the ball sections are certified. Floating-point peak searches are computed.
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
import json
import time
import signal
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent / "pylib"))
import mpmath as mp
import numpy as np
from scipy.special import digamma
from scipy.optimize import minimize_scalar
from scipy.integrate import quad
from flint import acb, acb_series, arb, ctx

signal.alarm(600)
mp.mp.dps = 60
ctx.dps = 65
ROOT = Path(__file__).parent
OUT = ROOT / "flow_results.json"
N = 2_000_000


def record(key, value):
    data = json.loads(OUT.read_text()) if OUT.exists() else {}
    data[key] = value
    OUT.write_text(json.dumps(data, indent=2) + "\n")
    print(key, json.dumps(value, indent=2), flush=True)


def primes_powers(n=N):
    sieve = np.ones(n + 1, dtype=bool)
    sieve[:2] = False
    for p in range(2, int(n**0.5) + 1):
        if sieve[p]:
            sieve[p*p::p] = False
    primes = np.flatnonzero(sieve)
    pairs = []
    for p in primes:
        q = int(p)
        while q <= n:
            pairs.append((q, int(p)))
            q *= int(p)
    pairs.sort()
    return primes, pairs


def fb(a, length):
    """Integral of exp(-a*x)*(1-x/L), with a=0 handled exactly."""
    if a == 0:
        return length / 2
    return 1/a - (1-(-a*length).exp())/(length*a*a)


def hb(z, delta, length):
    return fb(acb(delta)-acb(0,1)*z, length) + fb(acb(delta)+acb(0,1)*z, length)


def fm(a, length):
    if not a:
        return length / 2
    return 1/a + mp.expm1(-a*length)/(length*a*a)


def hm(z, delta, length):
    return fm(delta-1j*z, length)+fm(delta+1j*z, length)


def fn(a, length):
    return 1/a + np.expm1(-a*length)/(length*a*a)


def hn(z, delta, length):
    return fn(delta-1j*z, length)+fn(delta+1j*z, length)


def trigamma_np(z):
    z = np.asarray(z, dtype=complex)
    out = np.zeros_like(z)
    for _ in range(32):
        out += 1/(z*z)
        z = z + 1
    out += 1/z+1/(2*z*z)
    for k in range(1, 10):
        out += float(mp.bernoulli(2*k))/z**(2*k+1)
    return out


def arch_np(w, delta, length, gamma_shift=.5, conductor=1):
    a = delta + gamma_shift - 1j*np.asarray(w)
    q = a/2
    out = np.log(conductor/np.pi)+digamma(q).real+trigamma_np(q).real/(2*length)
    for n in range(3):
        b = a+2*n
        out -= 2/length*(np.exp(-b*length)/(b*b)).real
    return out


def arch_ball(w, delta, length, gamma_shift=None, conductor=1):
    gamma_shift = arb(1)/2 if gamma_shift is None else gamma_shift
    a = acb(delta+gamma_shift, -w)
    qseries = acb_series([a/2, 1], prec=3).lgamma()
    out = (arb(conductor)/arb.pi()).log()+qseries[1].real+qseries[2].real/length
    stop = 4
    for n in range(stop):
        b = a+2*n
        out -= 2/length*((-b*length).exp()/(b*b)).real
    first = delta+gamma_shift+2*stop
    tail = 2/length*(-first*length).exp()/(first*first*(1-(-2*length).exp()))
    out += arb(0, tail.abs_upper())
    return out


def ball_basics():
    primes, pp = primes_powers()
    length = arb(N).log()
    d = arb(1)/20
    w = arb("14.1347")
    prime = arb(0)
    logs = {}
    for n, p in pp:
        if p not in logs:
            logs[p] = arb(p).log()
        ln = arb(n).log()
        prime += 2*logs[p]*(-(arb(1)/2+d)*ln).exp()*(1-ln/length)*(w*ln).cos()
    pole = (hb(acb(w,arb(1)/2),d,length)+hb(acb(w,-arb(1)/2),d,length)).real
    arch = arch_ball(w,d,length)
    exact = pole+arch-prime
    vals = {"status":"computed", "arithmetic":"certified python-flint balls", "precision_decimal":ctx.dps,
            "prime_count":len(primes), "prime_power_count":len(pp), "L":str(length),
            "resolution_2pi_over_L":str(2*arb.pi()/length),
            "zeta_nzeros_100":str(arb(100).zeta_nzeros()),
            "zeta_nzeros_150":str(arb(150).zeta_nzeros()), "zeta_nzeros_160":str(arb(160).zeta_nzeros()),
            "sample_delta":"0.05", "sample_w":"14.1347", "sample_prime_sum":str(prime),
            "sample_pole_term":str(pole), "sample_arch_term_including_minus_log_pi":str(arch),
            "sample_O":str(exact), "reported_prime_value":"11.56902", "reported_zero_value":"11.56918",
            "H0":{ds:str(hb(acb(0),arb(ds),length).real) for ds in ["0.5","0.2","0.1","0.05"]}}
    record("ball_basics",vals)


def table():
    length = np.log(N)
    primes, pp = primes_powers()
    ns = np.array([n for n,p in pp],dtype=float)
    ps = np.array([p for n,p in pp],dtype=float)
    logs = np.log(ns)
    deltas = np.array([.5,.2,.1,.05])
    weights = (2*np.log(ps)*(1-logs/length))[:,None]*np.exp(-logs[:,None]*(deltas[None,:]+.5))
    grid = np.linspace(5,160,6201)
    values = np.empty((len(grid),4))
    for j in range(0,len(grid),32):
        ws = grid[j:j+32]
        prime = np.cos(ws[:,None]*logs[None,:]) @ weights
        for k,d in enumerate(deltas):
            pole = 2*hn(ws+.5j,d,length).real
            values[j:j+32,k] = pole+arch_np(ws,d,length)-prime[:,k]
    zballs = acb.zeta_zeros(1,59)
    zeros = np.array([float(z.imag.mid()) for z in zballs])
    rows=[]
    for k,d in enumerate(deltas):
        h0=float(hn(0j,d,length).real)
        ids=np.flatnonzero((values[1:-1,k]>values[:-2,k])&(values[1:-1,k]>values[2:,k]))+1
        ids=ids[values[ids,k]>.25*h0]
        def val(w):
            return float((2*hn(w+.5j,d,length).real+arch_np(w,d,length))-np.cos(w*logs) @ weights[:,k])
        positions=[]
        heights=[]
        for ix in ids:
            sol=minimize_scalar(lambda w:-val(w),bracket=(grid[ix-1],grid[ix],grid[ix+1]),method="brent",options={"xtol":1e-13,"maxiter":100})
            positions.append(float(sol.x))
            heights.append(float(-sol.fun/h0))
        positions=np.array(positions)
        errors=np.array([np.min(np.abs(positions-z)) for z in zeros if z<150])
        rows.append({"delta":float(d),"peaks_count":len(positions),"zeros_under_150":int(sum(zeros<150)),
                     "max_error_first52":float(errors.max()),"median_error_first52":float(np.median(errors)),
                     "normalized_height_min_all58":min(heights),"normalized_height_max_all58":max(heights),
                     "positions":positions.tolist(),"normalized_heights":heights,"errors_first52":errors.tolist()})
    record("table",{"status":"computed","arithmetic":"NumPy/SciPy float64, no interval extrema enclosure", "grid_step":.025,"rows":rows})


def dh_mp(s):
    k=(mp.sqrt(10-2*mp.sqrt(5))-2)/(mp.sqrt(5)-1)
    return mp.power(5,-s)*(mp.zeta(s,mp.mpf(1)/5)+k*mp.zeta(s,mp.mpf(2)/5)-k*mp.zeta(s,mp.mpf(3)/5)-mp.zeta(s,mp.mpf(4)/5))


def dh_ball(s):
    k=((10-2*arb(5).sqrt()).sqrt()-2)/(arb(5).sqrt()-1)
    coeff=[1,k,-k,-1]
    ans=s.zeta(arb(1)/5)
    for j in range(2,5):
        ans+=coeff[j-1]*s.zeta(arb(j)/5)
    return (-s*arb(5).log()).exp()*ans


def dh():
    root=mp.findroot(dh_mp,(mp.mpc("0.8085","85.6993"),mp.mpc("0.8086","85.6994")),solver="secant",tol=mp.mpf("1e-55"))
    cr=mp.nstr(root.real,55)
    ci=mp.nstr(root.imag,55)
    c=acb(cr,ci)
    r=arb("1e-20")
    box=c+acb(arb(0,r),arb(0,r))
    f0=dh_ball(c)
    deriv_center=dh_ball(acb_series([c,1],prec=2))[1]
    deriv_box=dh_ball(acb_series([box,1],prec=2))[1]
    q=(1-deriv_box/deriv_center).abs_upper()
    displacement=(f0/deriv_center).abs_upper()
    contraction=bool(q<1 and displacement+q*r<r)
    assert contraction
    assert box.real>arb(1)/2
    eta=box.real-arb(1)/2
    length=arb(N).log()
    ratio={}
    growth=[]
    for ds in ["0.5","0.2","0.1","0.05"]:
        d=arb(ds)
        R=(fb(acb(d-eta),length)+fb(acb(d+eta),length)).real/(2*fb(acb(d),length).real)
        ratio[ds] = str(R)
    for ls in [10,20,40,80,160]:
        d=arb("0.05")
        L=arb(ls)
        B=2*(fb(acb(d-eta),L)+fb(acb(d+eta),L)).real
        growth.append({"L":ls,"own_positive_height_pair_term":str(B)})
    negative_lobes=[]
    for shift in ["0","0.2","0.3","0.4","0.5","0.6","0.8","1"]:
        d=arb("0.05")
        B=(hb(acb(shift,eta),d,length)+hb(acb(shift,-eta),d,length)).real
        negative_lobes.append({"offset":shift,"pair_component":str(B),"normalized_component":str(B/hb(acb(0),d,length).real)})
    sb=acb("0.7","2.3")
    def completed(s):
        return (((s+1)/2)*(arb(5)/arb.pi()).log()).exp()*((s+1)/2).gamma()*dh_ball(s)
    fe=completed(sb)-completed(1-sb)
    vals={"status":"computed","arithmetic":"certified python-flint balls for contraction and displayed intervals",
          "mpmath_root":str(root), "center_real":cr,"center_imag":ci,"disk_radius":str(r),
          "D_center":str(f0),"Dprime_center":str(deriv_center),"Dprime_box":str(deriv_box),
          "contraction_q":str(q),"newton_displacement":str(displacement),
          "contraction_qr_plus_displacement":str(q*r+displacement),"strict_self_mapping":contraction,
          "root_containing_rectangle":str(box),"eta":str(eta),"height_ratio_R":ratio,
          "negative_lobes_at_manuscript_parameters":negative_lobes,
          "divergence_at_delta_005":growth,"functional_equation_sample_residual":str(fe),
          "no_zeros_Re_ge2_bound_zeta2_minus1":str(arb(2).zeta()-1)}
    record("davenport_heilbronn",vals)


def dh_flow():
    nmax=10000
    k=((10-2*arb(5).sqrt()).sqrt()-2)/(arb(5).sqrt()-1)
    aa=[arb(0),arb(1),k,-k,arb(-1)]
    coefficients=[arb(0)]*(nmax+1)
    for n in range(2,nmax+1):
        coefficients[n]=aa[n%5]*arb(n).log()
    for d in range(2,nmax//2+1):
        cd=coefficients[d]
        for m in range(2,nmax//d+1):
            coefficients[d*m]-=cd*aa[m%5]
    length=arb(nmax).log()
    delta=arb("0.05")
    ns=np.arange(2,nmax+1,dtype=float)
    lns=np.log(ns)
    coeff_np=np.array([float(c.mid()) for c in coefficients[2:]])
    weights=2*coeff_np*np.exp(-.55*lns)*(1-lns/float(length.mid()))
    grid=np.linspace(84.5,87,501)
    vals=arch_np(grid,.05,float(length.mid()),gamma_shift=1.5,conductor=5)-np.cos(grid[:,None]*lns[None,:]) @ weights
    minw=float(grid[np.argmin(vals)])
    maxw=float(grid[np.argmax(vals)])
    probes=[]
    for wstr in [repr(minw),repr(maxw),"85.69934848537759217192926770894"]:
        w=arb(wstr)
        total=arb(0)
        for n in range(2,nmax+1):
            logn=arb(n).log()
            total+=2*coefficients[n]*(-(arb(1)/2+delta)*logn).exp()*(1-logn/length)*(w*logn).cos()
        value=arch_ball(w,delta,length,gamma_shift=arb(3)/2,conductor=5)-total
        probes.append({"w":wstr,"O_D":str(value),"O_D_over_H0":str(value/hb(acb(0),delta,length).real)})
    record("dh_finite_flow",{"status":"computed","arithmetic":"certified balls for coefficients and probe values; floating-point grid chooses probes",
                             "N":nmax,"L":str(length),"delta":"0.05","lambda_D_6":str(coefficients[6]),
                             "lambda_D_6_formula":str((1+k*k)*arb(6).log()),"probes":probes})


def zeros():
    ctx.dps = 35
    zs = acb.zeta_zeros(1,400)
    d=arb("0.05")
    length=arb(N).log()
    w=arb("14.1347")
    total=arb(0)
    for z in zs:
        total+=(hb(acb(w-z.imag),d,length)+hb(acb(w+z.imag),d,length)).real
    T=float(zs[-1].imag.mid())
    dm=.05
    wm=14.1347
    lm=float(length.mid())
    # The smooth counting-density tail neglects S(t). Quadrature is not certified.
    def rho(t):
        return np.log(t/(2*np.pi))/(2*np.pi)
    def smooth(t):
        out=0.
        for r in [t-wm,t+wm]:
            den=dm*dm+r*r
            out+=2*dm/den-2/lm*(dm*dm-r*r)/(den*den)
        return rho(t)*out
    def cosine_coeff(t):
        rm,rp=t-wm,t+wm
        denm,denp=(dm*dm+rm*rm)**2,(dm*dm+rp*rp)**2
        return rho(t)*(np.cos(wm*lm)*((dm*dm-rm*rm)/denm+(dm*dm-rp*rp)/denp)+2*dm*np.sin(wm*lm)*(rm/denm-rp/denp))
    def sine_coeff(t):
        rm,rp=t-wm,t+wm
        denm,denp=(dm*dm+rm*rm)**2,(dm*dm+rp*rp)**2
        return rho(t)*(np.sin(wm*lm)*((dm*dm-rm*rm)/denm-(dm*dm-rp*rp)/denp)-2*dm*np.cos(wm*lm)*(rm/denm+rp/denp))
    main=quad(smooth,T,np.inf,epsabs=1e-13,limit=200)
    cosint=quad(cosine_coeff,T,np.inf,weight="cos",wvar=lm,epsabs=1e-13,limit=200)
    sinint=quad(sine_coeff,T,np.inf,weight="sin",wvar=lm,epsabs=1e-13,limit=200)
    tail=main[0]+2/lm*np.exp(-dm*lm)*(cosint[0]+sinint[0])
    record("zero_sum",{"status":"computed","first400_sum_ball":str(total),"T400":str(zs[-1].imag),
                       "smooth_integral_tail":tail,"sum_plus_smooth_tail":float(total.mid())+tail,
                       "quad_error_estimates":[main[1],cosint[1],sinint[1]],
                       "tail_is_certified":False})


if __name__=="__main__":
    start=time.time()
    mode=sys.argv[1] if len(sys.argv)>1 else "all"
    for name,func in [("ball",ball_basics),("dh",dh),("dhflow",dh_flow),("zeros",zeros),("table",table)]:
        if mode in (name,"all"):
            func()
    print("elapsed_seconds",round(time.time()-start,3),flush=True)
