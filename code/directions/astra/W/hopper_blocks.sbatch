#!/bin/bash
# Proposal only. This file has not been submitted.
#SBATCH --job-name=polya-block-fit
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=hopper_logs/blocks_%j.log
set -euo pipefail
cd /scratch/$USER/polya_spin_W
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
python3 -I -u fit_spins.py --families block --block-n 512,1024 --starts 12 --output-dir results_hopper_blocks
