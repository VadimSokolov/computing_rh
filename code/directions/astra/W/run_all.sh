#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
timeout 600 python3 -I fit_spins.py > numerical_run.log 2>&1
timeout 600 python3 -I fit_spins.py --families block --block-n 64,128,256 --starts 12 --output-dir results_extended > extended_run.log 2>&1
timeout 600 python3 -I validate_results.py --zeros > zero_run.log 2>&1
timeout 600 python3 -I validate_results.py > diagnostic_run.log 2>&1
timeout 30 python3 -I build_report.py
