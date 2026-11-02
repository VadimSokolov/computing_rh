"""Independent Davenport and Heilbronn control, with ball arithmetic.

mpmath locates a zero. Arb bounds on a disk then prove its existence and
uniqueness by Rouche's theorem. A small residual alone is not used as proof.
"""
from pathlib import Path
import platform
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent / "pylib"))
import flint
from flint import arb, acb, acb_series, ctx
import mpmath as mp
from referee_numeric import mobius_sieve


def main():
    ctx.dps = 80
    ctx.cap = 4
    mp.mp.dps = 75
    lines = []

    def emit(key, value):
        line = f"{key}: {value}"
        lines.append(line)
        print(line, flush=True)

    emit("python", platform.python_version())
    emit("python_flint", flint.__version__)
    emit("mpmath", mp.__version__)
    emit("precision_bits", ctx.prec)
    pi = arb.pi()
    kappa = ((10-2*arb(5).sqrt()).sqrt()-2)/(arb(5).sqrt()-1)
    km = (mp.sqrt(10-2*mp.sqrt(5))-2)/(mp.sqrt(5)-1)
    coeff = [arb(0), arb(1), kappa, -kappa, arb(-1)]
    emit("kappa", kappa)
    assert kappa > 0 and kappa < 1

    def Dmp(s):
        return mp.power(5, -s)*sum(a*mp.zeta(s, mp.mpf(r)/5)
                                   for r, a in [(1, 1), (2, km), (3, -km), (4, -1)])

    def D(s):
        if isinstance(s, acb_series):
            return (-s*acb(5).log()).exp()*sum(
                (coeff[r]*s.zeta(a=acb(r)/5, deflate=1) for r in range(1, 5)),
                acb_series([0], s.prec))
        return (-s*acb(5).log()).exp()*sum(
            (coeff[r]*s.zeta(acb(r)/5) for r in range(1, 5)), acb(0))

    def Lambda(s):
        return (((s+1)/2)*acb(5/pi).log()).exp()*((s+1)/2).gamma()*D(s)

    root = mp.findroot(Dmp, mp.mpc("0.81", "85.7"))
    emit("mpmath_root", mp.nstr(root, 72))
    c = acb(mp.nstr(root.real, 72), mp.nstr(root.imag, 72)).mid()
    r = arb(2)**(-80)
    box = c + acb(arb(0, r), arb(0, r))
    d = D(acb_series([c, 1], 2))[1].mid()
    derivative_box = D(acb_series([box, 1], 2))[1]
    residual = abs(D(c))
    variation = abs(derivative_box-d)
    lhs = residual + variation*r
    rhs = abs(d)*r
    assert lhs < rhs
    assert box.real > arb(1)/2
    emit("root_disk_center", c)
    emit("root_disk_radius", r)
    emit("D_center", D(c))
    emit("derivative_center", d)
    emit("derivative_box_deviation_bound", variation)
    emit("Rouche_left_bound", lhs)
    emit("Rouche_right_bound", rhs)
    emit("Rouche_ratio_upper", (lhs/rhs).upper())
    emit("root_real_enclosure", box.real)
    emit("root_imag_enclosure", box.imag)
    emit("complex_frequency_enclosure", (box-acb(arb(1)/2))/acb(0,1))
    emit("complex_spin_size_enclosure", pi/((box-acb(arb(1)/2))/acb(0,1)))

    for s in [acb("0.8", "2.3"), acb("0.3", 20), c]:
        a, b = Lambda(s), Lambda(1-s)
        delta = a-b
        assert delta.contains(0)
        emit(f"functional_equation_residual_s={s}", delta)
    emit("Euler_coefficient_a6", arb(1))
    emit("Euler_coefficient_a2_times_a3", -kappa*kappa)

    # Direct theta sum with an analytic tail enclosure. For t>=1/4 and
    # N=40 the function x*exp(-pi*t*x*x/5) decreases for x>=N.
    def H(t):
        N = 40
        b = pi*t/5
        value = sum((n*coeff[n % 5]*(-b*n*n).exp() for n in range(1, N+1)), arb(0))
        tail = (-b*N*N).exp()/(2*b)
        return value+arb(0, tail.upper())

    q = (-3*pi/5).exp()
    positive_margin = 2-1/(1-q)**2
    emit("kernel_positivity_margin", positive_margin)
    assert positive_margin > 0
    for t in [arb(1)/4, arb(1)/2, arb(1), arb(2), arb(4)]:
        delta = H(t)-t**(-arb(3)/2)*H(1/t)
        assert delta.contains(0)
        emit(f"theta_transform_residual_t={t}", delta)

    lhalf = Lambda(acb_series([acb(arb(1)/2), 1], 4))
    lone = Lambda(acb_series([1, 1], 4))
    emit("Lambda_half", lhalf[0])
    emit("Lambda_one", lone[0])
    emit("DH_level_one_mean", lone.log()[1])
    emit("DH_level_one_variance", 2*lone.log()[2])
    assert lone.log()[2].real > 0
    for u in [arb(0), arb(1)/2, arb(1), arb(2)]:
        kernel = 2*(3*u/2).exp()*H((2*u).exp())
        emit(f"DH_even_kernel_u={u}", kernel)
        fp = kernel*(u/2).exp()/lone[0].real
        fm = kernel*(-u/2).exp()/lone[0].real
        emit(f"DH_tilted_density_ratio_u={u}", fp/fm)
        assert (fp/fm-u.exp()).contains(0)

    # The same algebraic Mobius identity works for the off-line zero's
    # complex frequency. This is not a real random-variable construction.
    frequency = (c-acb(arb(1)/2))/acb(0, 1)
    z, theta = acb("0.7"), arb(1)/2
    exact = ((theta+z)**2+frequency**2)/(theta**2+frequency**2)
    mu = mobius_sieve(20000)
    value = acb(0)
    for k in range(1, 20001):
        if mu[k]:
            a = pi/(k*frequency)
            x, y = a*(theta+z), a*theta
            value += mu[k]*((x.sinh()/x)/(y.sinh()/y)).log()
    emit("DH_complex_per_zero_factor", exact)
    emit("DH_complex_Mobius_K20000", value.exp())
    emit("DH_complex_Mobius_K20000_error", value.exp()-exact)
    emit("DH_complex_per_zero_mean_formula", 2*theta/(theta**2+frequency**2))
    emit("DH_complex_per_zero_variance_formula",
         2*(frequency**2-theta**2)/(theta**2+frequency**2)**2)
    Path("results").mkdir(exist_ok=True)
    Path("results/referee_dh.txt").write_text("\n".join(lines)+"\n")


if __name__ == "__main__":
    main()
