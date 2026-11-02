import importlib
import sys
from pathlib import Path
sys.path.insert(0, '/private/tmp/d1lib')

print(sys.version)
for name in ['numpy', 'scipy', 'mpmath', 'flint']:
    try:
        module = importlib.import_module(name)
        print(name, getattr(module, '__version__', 'unknown'))
        if name == 'flint':
            print('arb methods', [s for s in dir(module.arb) if any(k in s for k in ['lerch', 'digamma', 'polygamma', 'zeta', 'rad', 'lower', 'upper'])])
            print('arb_mat methods', [s for s in dir(module.arb_mat) if any(k in s for k in ['chol', 'ldl', 'eig', 'inv'])])
    except ImportError as exc:
        print(name, str(exc))
for root in [Path('/private/tmp'), Path('/opt/homebrew/lib'), Path('.')]:
    if root.exists():
        print('candidates', root, [str(p) for p in root.iterdir() if any(s in p.name.lower() for s in ['env', 'python', 'flint', 'arb'])])
