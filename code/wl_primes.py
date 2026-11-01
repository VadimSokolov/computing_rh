# wl_primes.py: Work package 1, "prime dependence" (i): the Weil ladder of zeta with the prime sum truncated at Y < lambda^2 = X (Y = 1: no primes, only the archimedean and pole terms).
# Only the sum -2 sum_{n <= Y} Lambda(n) n^{-1/2} g(log n) is cut; the support [1/lambda, lambda] and the modes stay those of X. Built with the general routine of wl_dh.py (extends wl_weil.py).
# For each (X, N, Y): the least eigenvalue (mpmath at 3.5N + 60 digits and above), its sign certified in ball arithmetic (positive by a ball LDL^T, negative by a ball Rayleigh quotient of the
# eigenvector), the zeros of xi_hat below pi N/L, their distances to the zeta ordinates (err_log10 per ordinate) and the number matched to 1e-3.
# Usage (from data/ on Hopper): python3 wl_primes.py run --procs 40 --out wl_primes_scan.json     (scripts beside)
import argparse, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from multiprocessing import Pool
import wl_dh


def grid():
    tasks = []
    for X, Ns in ((13, (12, 20, 40)), (50, (20, 40)), (200, (20, 40)), (1000, (20, 40))):
        ys = sorted(set([1, 2, 3, 5, 7, 11, 13, 17, 23, 30, 50, 100, 300, 600, X] + [int(X ** t) for t in (0.5, 0.75)]))
        for N in Ns:
            yy = ys
            if N == 20 and X in (50, 200, 1000):
                yy = sorted(set(ys) | set(n for n, _ in wl_dh.W.prime_powers(X)))
            for Y in yy:
                if Y <= X:
                    tasks.append(dict(kind='zeta', parity='even', X=X, N=N, Y=Y, cert=True, roots=True))
    return tasks


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['run'])
    ap.add_argument('--procs', type=int, default=1)
    ap.add_argument('--out', default='wl_primes_scan.json')
    ap.add_argument('--max-cells', type=int, default=0)
    a = ap.parse_args()
    tasks = grid()
    tasks.sort(key=lambda t: -t['N'] * 10 - t['X'] / 1000)
    if a.max_cells:
        tasks = tasks[:a.max_cells]
    with Pool(a.procs) as p:
        res = list(p.imap_unordered(wl_dh._cell, tasks))
    json.dump(res, open(a.out, 'w'), indent=1)
