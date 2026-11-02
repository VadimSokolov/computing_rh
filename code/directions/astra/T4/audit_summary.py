"""Collect already computed evidence and reproducibility metadata."""
from pathlib import Path
import hashlib
import json
import platform
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent / "pylib"))
import flint
import mpmath
import numpy
import scipy

BASE=Path(__file__).resolve().parent
out=BASE/"results"
out.mkdir(exist_ok=True)
files=["flow_results.json","mobius_results.json","smallball_results.json",
       "ball_inversion_results.json","dh_smallball_results.json","control_results.json"]
scripts=["flow_audit.py","mobius_audit.py","smallball_audit.py",
         "ball_inversion.py","dh_smallball_control.py","control_checks.py","audit_summary.py"]
records=[]
for filename in files:
    value=json.loads((BASE/filename).read_text())
    records.append(filename+"\n"+json.dumps(value,indent=2)+"\n")
(out/"numbers.txt").write_text("\n".join(records))
manifest={"status":"computed","python":platform.python_version(),
          "versions":{"mpmath":mpmath.__version__,"python-flint":flint.__version__,
                      "numpy":numpy.__version__,"scipy":scipy.__version__},
          "sha256":{f:hashlib.sha256((BASE/f).read_bytes()).hexdigest()
                    for f in ["note.tex"]+scripts+files}}
(out/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
print("Wrote results/numbers.txt and results/manifest.json")
