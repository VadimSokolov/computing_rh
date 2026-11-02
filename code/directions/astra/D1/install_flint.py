"""Install the optional ball-arithmetic backend in this output folder."""
from pathlib import Path
import subprocess
import sys

target = Path(__file__).resolve().parent / '.deps'
subprocess.run([sys.executable, '-I', '-m', 'pip', 'install', '--disable-pip-version-check',
                '--no-input', '--only-binary=:all:', '--target', str(target),
                'python-flint'], check=True, timeout=25)
