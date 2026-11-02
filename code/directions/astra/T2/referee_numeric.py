"""Independent ball-arithmetic audit of note.tex. Run with Python 3.13.

PYTHONPATH=./pylib python3 referee_numeric.py
All finite numerical expressions use Arb through python-flint. The smooth
zero-density tail is explicitly an approximation to the omitted zero sum.
"""
from pathlib import Path
import hashlib
import platform
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent / "pylib"))
import flint
from flint import arb, acb, acb_series, ctx


def mobius_sieve(limit):
    mu = [1] * (limit + 1)
    primes = []
    composite = [False] * (limit + 1)
    mu[0] = 0
    for n in range(2, limit + 1):
        if not composite[n]:
            primes.append(n)
            mu[n] = -1
        for p in primes:
            if n * p > limit:
                break
            composite[n * p] = True
            if n % p == 0:
                mu[n * p] = 0
                break
            mu[n * p] = -mu[n]
    return mu


def main():
    ctx.dps = 60
    ctx.cap = 5
    lines = []

    def emit(key, value):
        line = f"{key}: {value}"
        lines.append(line)
        print(line, flush=True)

    emit("python", platform.python_version())
    emit("python_flint", flint.__version__)
    emit("precision_bits", ctx.prec)
    emit("note_sha256", hashlib.sha256(Path("note.tex").read_bytes()).hexdigest())
    pi = arb.pi()
    mu = mobius_sieve(20000)
    zeros = acb.zeta_zeros(1, 1000)
    gammas = [z.imag for z in zeros]
    assert all(z.real == arb(1)/2 for z in zeros)
    g = gammas[0]
    emit("gamma_1", g)
    emit("gamma_1000", gammas[-1])
    z2 = acb(arb(1)/2, 2*g).zeta()
    emit("zeta_at_half_plus_2_gamma_1_i", z2)
    assert not z2.contains(0)
    # Independently bracket a line zero using continuity and the functional
    # equation. This obstruction does not require knowing its zero index.
    center, radius = g.mid(), arb(2)**(-160)
    bracket = arb(center, radius)
    signs = []
    for label, t in [("lower", center-radius), ("upper", center+radius)]:
        s = acb(arb(1)/2, t)
        value = s*(s-1)/2 * (-s*acb(pi).log()/2).exp()*(s/2).gamma()*s.zeta()
        signs.append(value.real)
        emit(f"xi_at_zero_bracket_{label}", value)
    assert signs[0]*signs[1] < 0
    doubled_box_value = acb(arb(1)/2, 2*bracket).zeta()
    assert not doubled_box_value.contains(0)
    emit("line_zero_sign_bracket", bracket)
    emit("zeta_at_doubled_sign_bracket", doubled_box_value)

    def sinhc(w):
        return w.sinh() / w

    for z in [acb("0.7"), acb("0.2", 13)]:
        answer = ((z + arb(1)/2)**2 + g*g) / (arb(1)/4 + g*g)
        emit(f"per_zero_exact_z={z}", answer)
        log_product = acb(0)
        targets = {100, 2000, 20000} if z.imag == 0 else {20000}
        for k in range(1, 20001):
            if mu[k]:
                a = pi / (k*g)
                term = sinhc(a*(z+arb(1)/2))/sinhc(acb(a/2))
                log_product += mu[k]*term.log()
            if k in targets:
                value = log_product.exp()
                emit(f"per_zero_K={k}_z={z}", value)
                emit(f"per_zero_error_K={k}_z={z}", value-answer)
                # For k>K all Euler factors satisfy |w| <= q < 1.
                # Bound the absolute log of the omitted product using zeta(2).
                q = max(abs(z+arb(1)/2).upper(), arb(1)/2)**2 / ((k+1)*g)**2
                bound = ((abs(z+arb(1)/2)**2+arb(1)/4)
                         * arb(2).zeta()/(g*g*k*(1-q)))
                emit(f"absolute_log_tail_bound_K={k}_z={z}", bound)

    # Analytic continuation at s=1 using deflated zeta, no finite difference.
    x = acb_series([0, 1], 5)
    s = 1+x
    xi = s/2 * (1+x*s.zeta(deflate=1)) * (-s*acb(pi).log()/2).exp() * (s/2).gamma()
    logxi = xi.log()
    direct_mean = logxi[1].real
    direct_var = 2*logxi[2].real
    emit("direct_mean", direct_mean)
    emit("direct_variance", direct_var)

    # The finite 1000 by 3000 spin sums are evaluated through a convergent
    # partial-fraction series, with explicit absolute remainder bounds.
    K, J = 3000, 16
    theta = arb(1)/2
    powers = [sum((gg**(-2*j) for gg in gammas), arb(0)) for j in range(1, J+2)]
    dirichlet = [sum((arb(mu[k])/arb(k)**(2*j) for k in range(1, K+1)), arb(0))
                 for j in range(1, J+1)]
    sm, sv = arb(0), arb(0)
    for j in range(1, J+1):
        term = 2 * (-1)**(j-1) * arb(2*j).zeta()*powers[j-1]*dirichlet[j-1]
        sm += term*theta**(2*j-1)
        sv += term*(2*j-1)*theta**(2*j-2)
    q = (theta/g)**2
    rm = 2*theta**(2*J+1)*arb(2*J+2).zeta()**2*powers[J]/(1-q)
    rv = (2*theta**(2*J)*arb(2*J+2).zeta()**2*powers[J]
          *((2*J+1)/(1-q)+2*q/(1-q)**2))
    sm += arb(0, rm.upper())
    sv += arb(0, rv.upper())
    emit("spin_mean_1000_K3000", sm)
    emit("spin_variance_1000_K3000", sv)
    emit("spin_series_mean_remainder_bound", rm)
    emit("spin_series_variance_remainder_bound", rv)
    em = sum((1/(arb(1)/4+gg*gg) for gg in gammas), arb(0))
    ev = sum((2*(gg*gg-arb(1)/4)/(gg*gg+arb(1)/4)**2 for gg in gammas), arb(0))
    emit("exact_factor_mean_1000", em)
    emit("exact_factor_variance_1000", ev)
    emit("spin_mean_K_error_1000", sm-em)
    emit("spin_variance_K_error_1000", sv-ev)

    # Evaluate smooth-density tail integrals by an alternating power series.
    # These intervals enclose the integrals, not the true omitted-zero tails.
    gnext = acb.zeta_zero(1001).imag
    emit("gamma_1001", gnext)
    for label, T in [("at_gamma1000", gammas[-1]),
                     ("at_midpoint1000_1001", (gammas[-1]+gnext)/2)]:
        L = (T/(2*pi)).log()
        tm, tv = arb(0), arb(0)
        for j in range(9):
            p = 2*j+1
            base = (-arb(1)/4)**j / T**p * (L/p+arb(1)/p**2)
            if j < 8:
                tm += base/(2*pi)
                tv += p*base/pi
            else:
                tm += arb(0, abs(base/(2*pi)).upper())
                tv += arb(0, abs(p*base/pi).upper())
        emit(f"smooth_cutoff_{label}", T)
        emit(f"smooth_mean_tail_integral_{label}", tm)
        emit(f"smooth_variance_tail_integral_{label}", tv)
        emit(f"spin_plus_smooth_mean_{label}", sm+tm)
        emit(f"spin_plus_smooth_variance_{label}", sv+tv)
        emit(f"spin_plus_smooth_mean_error_{label}", sm+tm-direct_mean)
        emit(f"spin_plus_smooth_variance_error_{label}", sv+tv-direct_var)

    # Conditional-on-RH distances for the untilted sums, via an explicit
    # coupling. Q = sum gamma^-2 is computed unconditionally from xi at 1/2;
    # its interpretation as a real zero sum uses RH.
    s = acb_series([acb(arb(1)/2), 1], 5)
    xi_half = s*(s-1)/2 * (-s*acb(pi).log()/2).exp() * (s/2).gamma()*s.zeta()
    Q = (xi_half.log()[2]).real
    emit("Q_half_log_coefficient", Q)
    A = arb(2).zeta()/arb(4).zeta()
    B = 1/arb(2).zeta()
    Cminus = (A-B)/2
    Cextra = (A+B)/2-1
    factor = pi*pi/3*Q
    vm = factor*Cminus
    vr = factor*Cextra
    emit("mobius_negative_weight", Cminus)
    emit("mobius_positive_k_gt_1_weight", Cextra)
    emit("negative_variance_RH", vm)
    emit("positive_extra_variance_RH", vr)
    emit("k2_fraction_of_negative_variance", (arb(1)/4)/Cminus)
    emit("W2_upper_bound_RH", vm.sqrt()+vr.sqrt())
    emit("W2_lower_bound_RH", (factor.sqrt()-(2*Q).sqrt()))

    # Small exact arithmetic counterexamples for the frequency-count claims.
    for ordinates, point in [([2, 3], 6), ([2, 2], 2), ([2, 4], 8)]:
        c = sum(mu[point//g0] for g0 in ordinates if point % g0 == 0)
        emit(f"aggregate_c_ordinates={ordinates}_x={point}", c)
    Path("results").mkdir(exist_ok=True)
    Path("results/referee_numeric.txt").write_text("\n".join(lines)+"\n")


if __name__ == "__main__":
    main()
