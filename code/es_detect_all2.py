# es_detect_all2.py: runs es_detect.py for the five noise levels and two seeds in parallel (Hopper driver). Output: one JSON line per run on stdout.
import subprocess, sys, os
from concurrent.futures import ThreadPoolExecutor
here = os.path.dirname(os.path.abspath(__file__))
def go(a):
    t, s = a
    r = subprocess.run([sys.executable, os.path.join(here, 'es_detect.py'), '--tau', t, '--nnull', '1200', '--npow', '150', '--seed', s, '--om0', '2', '--Om', '21.8'], capture_output=True, text=True, cwd=here, env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1'))
    return r.stdout.strip() + r.stderr.strip()[-300:]
with ThreadPoolExecutor(8) as ex:
    for o in ex.map(go, [(t, s) for s in ('1', '2') for t in ('1e-10', '1e-8', '1e-6', '1e-4', '1e-2')]): print(o, flush=True)
