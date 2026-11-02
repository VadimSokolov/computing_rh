"""Independent first-rung pilot. All values are floating point, not certificates."""
import json
import math
from pathlib import Path
import numpy as np
from scipy import special, linalg
from numpy.polynomial.legendre import leggauss


def lower_matrix(L, T, nmax, order=12, panel_width=0.25):
    a = math.log(2)
    c = math.sqrt(2) * a
    beta = math.log(T / (2 * math.pi)) - 1 / (3 * T) - 0.5 / T**2 - c
    z, w = leggauss(order)
    edges = np.linspace(0, T, math.ceil(T / panel_width) + 1)
    mids, halves = (edges[1:] + edges[:-1]) / 2, np.diff(edges) / 2
    t = (mids[:, None] + halves[:, None] * z).ravel()
    weights = (halves[:, None] * w).ravel()
    n = np.arange(nmax)
    b = (np.sqrt(2 * L * (2 * n + 1))[:, None]
         * special.spherical_jn(n[:, None], L * t[None, :])
         * (-1.)**(n[:, None] // 2))
    m = special.digamma(0.25 + 0.5j * t).real - math.log(math.pi)
    pole = np.sqrt(2 * L * (2 * n + 1)) * special.spherical_in(n, L / 2)
    out = {}
    for parity in (0, 1):
        idx = np.arange(parity, nmax, 2)
        v = b[idx]
        pp = pole[idx]
        prime = (v * (weights * c * np.cos(a * t) / math.pi)) @ v.T
        arch = beta * np.eye(len(idx)) + (v * (weights * (m-beta) / math.pi)) @ v.T
        arch += (2 if parity == 0 else -2) * np.outer(pp, pp)
        M = arch - prime
        eig, vec = linalg.eigh(M)
        aeig, avec = linalg.eigh(arch)
        # The prime matrix here is band truncated, for diagnostics only.
        principal = avec[:, 1:]
        comp = linalg.eigvalsh(principal.T @ M @ principal)
        out[str(parity)] = {
            'lambda_min_R': float(eig[0]),
            'lambda_2_R': float(eig[1]),
            'arch_first': [float(x) for x in aeig[:3]],
            'R_on_arch_ground_complement': float(comp[0]),
            'prime_Rayleigh_at_bottom': float(vec[:, 0] @ prime @ vec[:, 0]),
        }
    # Uniform transform tail beyond order nmax via Rodrigues' formula.
    r = (L*T)**2 / ((2*nmax+1)*(2*nmax+3))
    loga2 = math.log(2*L*(2*nmax+1)) + 2*nmax*math.log(L*T)
    loga2 -= 2*(special.gammaln(2*nmax+2)-nmax*math.log(2)-special.gammaln(nmax+1))
    S = math.exp(loga2)/(1-r)
    # Safe real-band bound: digamma real part increases with t.
    W = max(abs(special.digamma(0.25)-math.log(math.pi)-beta),
            abs(special.digamma(0.25+0.5j*T).real-math.log(math.pi)-beta))+c
    eC = W*T/math.pi*(2*math.sqrt(2*L*S)+S)
    rp = (L/2)**2/((2*nmax+1)*(2*nmax+3))
    logp2 = math.log(2*L*(2*nmax+1))+L+2*nmax*math.log(L/2)
    logp2 -= 2*(special.gammaln(2*nmax+2)-nmax*math.log(2)-special.gammaln(nmax+1))
    Sp = math.exp(logp2)/(1-rp)
    ep = 2*(2*math.sqrt(2*math.sinh(L)*Sp)+Sp)
    return {'L': L, 'T': T, 'nmax': nmax, 'quadrature_order': order,
            'beta': beta, 'W': float(W), 'tail_S': S,
            'operator_truncation_bound': eC+ep, 'parities': out}


if __name__ == '__main__':
    L = math.log(3)/2
    rows = [lower_matrix(L, T, n) for T, n in [(25, 64), (40, 80), (60, 64), (60, 96), (100, 128)]]
    rows.append(lower_matrix(L, 40, 80, order=20))
    output = {'status': 'computed, not certified', 'endpoint': 'L=log(3)/2', 'runs': rows}
    Path('pilot_results.json').write_text(json.dumps(output, indent=2)+'\n')
    for row in rows:
        print(json.dumps(row))
