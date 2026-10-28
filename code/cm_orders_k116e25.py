# cm_orders_k116e25.py: complete monotonicity of the heat trace at every order on (0, 1.70], for the subsection "The prime side of the ladder" (ch/ch12.tex).
# Part (b) of code/cm_orders.py with K = 1.16 x 10^25 (the orders certified by Theorem thm:ch12:window) in place of K = 10^17, at the time T1 = 1.70, and with the
# prime powers n <= 10^9 in ball arithmetic (cm_orders.py stops at 2^24, which leaves a tail bound of about 1.1 at t = 1.7, larger than the margin). Same inequality,
# 2 c_k(t)/M_k(t) >= E_k(t) - S_C(t) - P_k(t) for every k > K and t in [T0, T1]; same lemma constants (part (b) of the lemma, eps = 10^-6); the part (a) of cm_orders.py (t <= T0) is not repeated.
# Every inequality is evaluated in ball arithmetic (python-flint, Arb) and counts only if it holds for the whole ball.
# Run: python3 cm_orders_k116e25.py [T1]; writes cm_orders_k116e25.json.
import json
import math
import sys
import time

import mpmath as mp
import numpy as np
from flint import arb, acb, ctx

ctx.prec = 128
pi, half, E = arb.pi(), arb(1) / 2, arb(1).exp()
K = 116 * 10**23
T0 = arb(19) / 1000
out, holds = {'certificate': {}, 'float': {}}, []
T00 = time.time()


def s(x, n=12):
    return x.str(n, radius=False) if isinstance(x, arb) else x


def check(name, cond, **vals):
    holds.append(bool(cond))
    out['certificate'][name] = dict(holds=bool(cond), **{k: s(v) for k, v in vals.items()})
    print(name, bool(cond), round(time.time() - T00, 1), flush=True)


def m(r):  # m(r) = Re psi(1/4 + i r/2) - log pi, for a real ball r
    return acb(arb(1) / 4, r / 2).digamma().real - pi.log()


C1, C2T = arb('2.62'), arb('2.08')
c2 = (pi * (2 * E + half)).sqrt().sqrt()
check('tail_constant', c2 < C2T, c2=c2)

NB = 10**9
EPSJ = arb(1) / 10**6
K1 = arb(K + 1)
C_large = (((K1 + arb('0.37')) / ((1 - EPSJ) * K1 + arb(1) / 4)).sqrt()
           + 2 * (K1 + arb('0.37')).sqrt() * (-EPSJ * K1).exp()).upper()
X_large = 2 * (EPSJ * K1).sqrt() / E
x_cut = arb(NB).log() / (2 * T0.sqrt())
check('b_lemma_monotone', (arb(1) / 4 < arb('0.37') * (1 - EPSJ)) and (K1 > 1 / (2 * EPSJ)))
check('b_lemma_constant', C_large < 1 + 6 * arb(10) ** -7, C=C_large)
check('b_lemma_range', x_cut < X_large, x_at_cutoff=x_cut, X=X_large)


def prime_power_blocks(N):  # arrays (n, p) of the prime powers n = p^j <= N, in blocks
    sieve = np.ones(N + 1, bool)
    sieve[:2] = False
    for p in range(2, int(N**0.5) + 1):
        if sieve[p]:
            sieve[p * p::p] = False
    ps = np.nonzero(sieve)[0]
    del sieve
    for i in range(0, len(ps), 10**6):
        yield ps[i:i + 10**6].tolist()


def SC_up(t, N):  # an upper bound for S_C(t): the constant C for the prime powers n <= N, the tail n > N bounded as in cm_orders.py
    lt, total, count = 2 * t.sqrt(), arb(0), 0
    for block in prime_power_blocks(N):
        for p in block:
            lp = arb(p).log()
            q_ = p
            while q_ <= N:
                ln = arb(q_).log()
                total += 2 * lp / arb(q_).sqrt() * (-(ln * ln) / (8 * t)).exp()
                count += 1
                q_ *= p
        print('  blocks done, last prime', block[-1], round(time.time() - T00, 1), flush=True)
    out['prime_powers'] = count
    q = arb(N).log() / (8 * t)
    assert N >= math.exp(10) and q > 1 and arb(N).log() / lt >= 1
    tail = (2 * (C1 + C2T / lt.sqrt()) * arb(N) ** (1 - q) / (q - 1)).upper()
    return (C_large * total + tail).upper(), tail


mp.mp.dps = 30


def EK_f(t):  # E_{K+1}(t) = (psi(K + 3/2) - log t)/2 - log 2 pi, up to a relative correction of order t/K
    return float((mp.digamma(K + 1.5) - mp.log(t)) / 2 - mp.log(2 * mp.pi))


alpha, eps = arb(K) + arb(3) / 2, arb(1) / 8192
a = alpha * (1 - eps) ** 2
logp = alpha * (2 * (1 - eps).log() + 1 - (1 - eps) ** 2)
check('b_chernoff', logp < -1000, log_bound=logp)
m0 = m(arb(0))
Pb = arb(-1000).exp().upper()
T1 = arb(sys.argv[1] if len(sys.argv) > 1 else '1.70')
mr1 = m((a / T1).sqrt())
eK = (mr1 - (mr1 - m0) * arb(-1000).exp()).lower()
Sb, tailb = SC_up(T1, NB)
check('b_T1_between_T0_and_2', (T0 < T1) and (T1 <= 2), T1=T1)
logP = (4 * pi).log() - (K + 1) * arb(4).log() + T1 / 4 + (K + arb(3) / 2) * T1.log() - alpha.lgamma()
check('b_polar_small', logP < -1000, log_P=logP)
check('b_tail_below_1e-3', tailb < arb(10) ** -3, S_tail=tailb)
check('b_margin', eK - Sb - Pb > 0, T1=T1, E_K1_lower=eK, S_C_upper=Sb, S_tail=tailb, P_upper=Pb, margin=eK - Sb - Pb,
      E_K1_float=EK_f(float(T1.mid())))
out['K'] = str(K)
out['T1'] = s(T1)
out['all_hold'] = all(holds)
out['seconds'] = round(time.time() - T00, 1)
print('certificate', all(holds), 'T1', s(T1), 'seconds', out['seconds'], flush=True)
json.dump(out, open('cm_orders_k116e25.json', 'w'), indent=1)
