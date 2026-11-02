import importlib
import json

out = {}
for name in ['numpy', 'scipy', 'mpmath', 'flint']:
    try:
        mod = importlib.import_module(name)
        out[name] = getattr(mod, '__version__', 'available')
    except ImportError as exc:
        out[name] = str(exc)
print(json.dumps(out, indent=2))
