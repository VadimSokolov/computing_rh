# zh_zeros.py: the ordinates of the first 8000 zeros of zeta (heights up to 8000th zero, about 10^4) as 50 digit decimal strings, one per line, in data/zeros_8000.txt, for code/zh_detect.py (the detection experiment of Chapter 18), whose sensitivity runs use 4000 and 8000 zeros; data/zeros_1700.txt stops at height 2197.26. Added for the book in October 2026, written from code/zeros_hp.py (same method).
# Method. python-flint (Arb) acb.zeta_zeros isolates and refines the zeros rigorously at 192 bits, in chunks of 1000. Checks: every ball has real part exactly 1/2 and imaginary radius below 1e-50; the ordinates increase; the strings reproduce the balls to a relative 1e-49; the first zeros agree with data/zeros_1700.txt (70 digits) to 1e-45 in absolute value; and the rigorous counting function arb.zeta_nzeros just above the last ordinate encloses the number of zeros written (when the enclosure is not a single integer the script says so). Nothing is rounded before the checks.
# Output: data/zeros_8000.txt (50 significant digits per line) and a log on standard output. Usage: python3 zh_zeros.py [--count 8000] [--prec 192] [--digits 50] [--out FILE] [--check FILE] [--quick]; --quick computes 60 zeros, checks them against zeros_1700.txt and writes to /tmp unless --out is given (about 1 s); the full case takes about one minute and is meant for the Hopper cluster (misc/reviews/nick-zhu/z2-run/zeros.slurm).
import argparse
import time
from pathlib import Path

from flint import acb, arb, ctx

HERE = Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--count', type=int, default=8000)
    ap.add_argument('--prec', type=int, default=192)
    ap.add_argument('--digits', type=int, default=50)
    ap.add_argument('--chunk', type=int, default=1000)
    ap.add_argument('--out')
    ap.add_argument('--check')
    ap.add_argument('--quick', action='store_true')
    args = ap.parse_args()
    count = 60 if args.quick else args.count
    out = args.out or ('/tmp/zeros_quick.txt' if args.quick else str(HERE.parent / 'data' / 'zeros_8000.txt'))
    ctx.prec = args.prec
    half = arb(1) / 2
    t0 = time.time()
    S, last = [], None
    for start in range(1, count + 1, args.chunk):
        n = min(args.chunk, count - start + 1)
        Z = acb.zeta_zeros(start, n)
        assert len(Z) == n
        for z in Z:
            assert z.real == half, 'a real part is not exactly 1/2'
            assert float(z.imag.rad()) < 1e-50, z.imag.rad()
            if last is not None:
                assert z.imag > last, 'ordinates not increasing'
            last = z.imag
            s = z.imag.str(args.digits, radius=False)
            assert abs(arb(s) - z.imag) < arb('1e-49') * z.imag, s
            S.append(s)
        print('zeros %d to %d done, %.1f s, gamma_%d = %s' % (start, start + n - 1, time.time() - t0, start + n - 1, S[-1]), flush=True)
    with open(out, 'w') as fh:
        fh.write('\n'.join(S) + '\n')
    print('wrote %d zeros to %s' % (len(S), out))
    chk = Path(args.check) if args.check else None
    if chk is None:
        for c in (HERE.parent / 'data' / 'zeros_1700.txt', HERE / 'zeros_1700.txt', Path('zeros_1700.txt')):
            if c.is_file():
                chk = c
                break
    if chk is not None:
        ref = [s.strip() for s in open(chk) if s.strip()]
        m = min(len(ref), len(S))
        d = max(abs(arb(ref[i]) - arb(S[i])) for i in range(m))
        print('compared with %s: %d zeros, max |difference| = %s' % (chk, m, d.str(5, radius=False)))
        assert d < arb('1e-45'), 'disagreement with the reference zeros'
    tt = arb(S[-1]) + arb('0.001')
    try:
        nz = arb.zeta_nzeros(tt)
        print('arb.zeta_nzeros(gamma_%d + 0.001) = %s (written: %d)' % (len(S), nz.str(10), len(S)))
    except Exception as e:
        print('zeta_nzeros not available or failed: %s' % e)
    print('done in %.1f s' % (time.time() - t0))


if __name__ == '__main__':
    main()
