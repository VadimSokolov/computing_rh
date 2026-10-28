# Note W (Chapter 13, Section fm:sec:bessel, paragraph "A small distance does not pin the zeros"), computation N5:
# the mean and variance of the clock T_q of the second row of Table fm:tab:bessel (Theorem fm:thm:dual, with b = 0):
#   E T_q = m_2/2,  Var T_q = -kappa_4/12,  kappa_4 = m_4 - 3 m_2^2,  m_k = int_R t^k q(t) dt,
# for q = m(t) exp(-a cosh l t)/C with the exact decimals of row 2 in wf_fit.json (book units). Each m_k = (2/C) int_0^inf
# t^k m exp(-a cosh l t) dt is enclosed by Arb's quadrature on [0, 5] plus the explicit tail bound
# int_5^inf t^k m exp(-a cosh l t) <= g(5)/h(5), g = exp(Lam t + k t - (a/2) e^{l t}) (wf_tvcert.Member.tail with k), and
# C is enclosed by the Bessel formula (wf_tvcert.Member.C_bessel), precision 192 bits. Check of the identities and of the
# quadrature: the same formulas for the Polya law p (m_k = (2/Xi(0)) int_0^inf t^k Phi, Phi summed to n = 8 with the
# remainder of fm_tvcert.Phi_ser, tail by fm_tvcert.tail_Phi) must give E T = 0.0231049931 and Var T = 3.71725993e-5.
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
# Usage: python3 wf_clock.py [--quick]  reads wf_fit.json; writes wf_clock.json
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fm_tvcert as T
import wf_tvcert as W
from flint import arb, acb, ctx

PREC = 192
ctx.prec = PREC
fit = json.load(open('wf_fit.json'))
mem = W.Member.from_row('row2', fit['row2']['exact'])
assert mem.mu == 0


def s6(x):
    return x.str(6, radius=False)


def moments_q():
    C = mem.C_bessel()
    Cq = mem.C_quad()
    assert C.overlaps(Cq)
    m = {}
    for k in (0, 2, 4):
        I = acb.integral(lambda z, _: z**k*mem.f(z), 0, 5, rel_tol=2**-(PREC - 30))
        m[k] = 2*(I.real + arb(0, mem.tail(5, k).upper()))/C
    return m, C


def moments_p():
    T.TB = (4*(2*arb.pi()**2*9**4 + 3*arb.pi()*81)*(2*arb.pi()*81 + arb('6.5'))**3*(-81*arb.pi()).exp()).upper()
    X0 = T.xi_half()
    def Phi(z):
        s = acb(0); PI = arb.pi()
        for n in range(1, 9):
            s += ((z*arb('4.5')).exp()*(2*PI**2*n**4) - (z*arb('2.5')).exp()*(3*PI*n**2))*((z*2).exp()*(-PI*n**2)).exp()*2
        return s + arb(0, T.TB)
    m = {}
    for k in (0, 2, 4):
        I = acb.integral(lambda z, _: z**k*Phi(z), 0, 3, rel_tol=2**-(PREC - 30))
        m[k] = 2*(I.real + arb(0, T.tail_Phi(k).upper()))/X0
    return m, X0


def clock(m):
    ET = m[2]/2
    k4 = m[4] - 3*m[2]**2
    return ET, -k4/12, k4


if __name__ == '__main__':
    t0 = time.time()
    mq, C = moments_q()
    ETq, VTq, k4q = clock(mq)
    mp_, X0 = moments_p()
    ETp, VTp, k4p = clock(mp_)
    out = dict(prec=PREC, row2=fit['row2']['exact'], C=C.str(20),
               q=dict(m0=mq[0].str(15), m2=mq[2].str(15), m4=mq[4].str(15), kappa4=k4q.str(12), ET=ETq.str(15), VarT=VTq.str(12),
                      ET_6=s6(ETq), VarT_6=s6(VTq)),
               p_check=dict(m0=mp_[0].str(15), m2=mp_[2].str(15), ET=ETp.str(15), VarT=VTp.str(12), ET_6=s6(ETp), VarT_6=s6(VTp),
                            book_ET='0.0231049931', book_VarT='3.71725993e-5'),
               secs=round(time.time() - t0, 1))
    print(json.dumps(out, indent=1))
    json.dump(out, open('wf_clock.json', 'w'), indent=1)
