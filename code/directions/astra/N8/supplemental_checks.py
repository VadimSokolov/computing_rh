"""Independent prime formula and normalization checks, with no source-code imports.

Run: timeout 30 python3 -I supplemental_checks.py
"""
import json
from pathlib import Path
import mpmath as mp

mp.mp.dps = 50
OUT = Path(__file__).resolve().parent


def xi(s):
    return s*(s-1)*mp.power(mp.pi,-s/2)*mp.gamma(s/2)*mp.zeta(s)/2


def mangoldt_terms(limit):
    sieve = bytearray(b'\x01')*(limit+1)
    sieve[:2] = b'\x00\x00'
    terms = []
    for p in range(2,limit+1):
        if sieve[p]:
            if p*p <= limit:
                sieve[p*p:limit+1:p] = b'\x00'*((limit-p*p)//p+1)
            power = p
            while power <= limit:
                terms.append((power,mp.log(p)))
                power *= p
    return sorted(terms)


terms = mangoldt_terms(int(mp.exp(12)))


def psi_primes(t):
    t = mp.mpf(t)
    finite = sum(weight/mp.sqrt(n)*(t-mp.log(n)) for n,weight in terms if mp.log(n) <= t)
    smooth = (4*(mp.exp(t/2)+mp.exp(-t/2)-2)
              +t/2*(mp.digamma(mp.mpf(1)/4)-mp.log(mp.pi))
              +(mp.pi**2+8*mp.catalan-mp.exp(-t/2)*mp.lerchphi(mp.exp(-2*t),2,mp.mpf(1)/4))/4)
    return smooth-finite


points = ['0.1','0.464002','1','2','8.039316023094285','12']
data = {"precision_dps":mp.mp.dps, "status":"computed, not certified",
        "psi_from_primes":{t:str(psi_primes(t)) for t in points},
        "number_prime_power_terms_up_to_exp12":len(terms),
        "bridge_transform_lambda1":str(mp.sqrt(2)/mp.sinh(mp.sqrt(2))),
        "squared_bridge_transform_lambda1":str((mp.sqrt(2)/mp.sinh(mp.sqrt(2)))**2),
        "multiplicity_two_at_gamma1":{
            "scalar_energy_variance_coefficient":4,
            "clock_variance_coefficient":2,
            "unit":"gamma^(-4)"}}
if (OUT/'independent_zeta_zeros.json').exists():
    zeros = [mp.mpf(s) for s in json.loads((OUT/'independent_zeta_zeros.json').read_text())]
    data['zero_count_used']=len(zeros)
    data['psi_partial_at_same_points']={t:str(sum(2*(1-mp.cos(g*mp.mpf(t)))/g**2 for g in zeros)) for t in points}
    if len(zeros) == 400:
        candidates = [(n,g,mp.siegeltheta(g)/mp.pi) for n,g in enumerate(zeros,1) if 30<g<680]
        low = min(candidates,key=lambda item:item[0]-2-item[2])
        high = max(candidates,key=lambda item:item[0]-1-item[2])
        data['phase_witnesses'] = {
            'minimum_left_limit':{'index':low[0],'height':str(low[1]),'S':str(low[0]-2-low[2])},
            'maximum_right_limit':{'index':high[0],'height':str(high[1]),'S':str(high[0]-1-high[2])},
            'S_oscillation':str(high[0]-1-high[2]-(low[0]-2-low[2])),
            'rank_one_all_phases_margin_over_two':str(high[0]-1-high[2]-(low[0]-2-low[2])-2)}
        full_mean=mp.diff(xi,mp.mpf('.5'),2)/xi(mp.mpf('.5'))
        data['RH_uniform_screw_tail_bound']=str(2*(full_mean-sum(2/g**2 for g in zeros)))
data['laplace_identity_values']={}
for z in [1,2,3]:
    left = mp.diff(xi,mp.mpf('.5')+z)/xi(mp.mpf('.5')+z)/z**2
    data['laplace_identity_values'][str(z)]=str(left)
(OUT/'supplemental_numbers.json').write_text(json.dumps(data,indent=2)+'\n')
print(json.dumps(data,indent=2))
