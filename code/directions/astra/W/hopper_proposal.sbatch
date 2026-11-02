#!/bin/bash
# Proposal only. This file has not been submitted.
#SBATCH --job-name=polya-spin-fit
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --array=0-95%8
#SBATCH --output=hopper_logs/fit_%A_%a.log
set -euo pipefail
cd /scratch/$USER/polya_spin_W
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
python3 -I -u hopper_fit.py --task "$SLURM_ARRAY_TASK_ID"
