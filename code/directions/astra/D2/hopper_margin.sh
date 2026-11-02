#!/usr/bin/env bash
# Proposal only. Arguments: staged wl_weil.py path, output directory.
#SBATCH --job-name=d2-margin
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=d2-margin-%j.log
set -euo pipefail
if [ "$#" -ne 2 ]; then
  echo "Usage: sbatch hopper_margin.sh /path/to/wl_weil.py /path/to/output"
  exit 2
fi
D2_WEIL_SCRIPT="$1"
D2_RUN_DIR="$2"
mkdir -p "$D2_RUN_DIR"
python3 -I "$D2_WEIL_SCRIPT" cert \
  --tasks '64:32,64:64,128:128,256:128,256:256,512:512' \
  --procs 1 \
  --out "$D2_RUN_DIR/d2-boundary-cells.json"
