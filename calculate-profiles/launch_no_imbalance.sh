#!/bin/bash
# Run this to submit: bash launch_no_imbalance.sh
# T = 0, 0.5, 0.9 (fixed); zeta = 0; int = -2:0.2:3

# --- Define ranges here ---
INT_START=-2.0; INT_STOP=3.0; INT_STEP=0.2
N_T=3

# --- Auto-calculate counts ---
N_INT=$(awk "BEGIN { print int(($INT_STOP - $INT_START) / $INT_STEP + 0.5) + 1 }")
N_JOBS=$(( N_T * N_INT ))

OUTDIR=out_no_imbalance

echo "Submitting $N_T x $N_INT = $N_JOBS jobs"

mkdir -p $OUTDIR

sbatch \
  --array=0-$(( N_JOBS - 1 )) \
  --output=/dev/null \
  --export=ALL,OUTDIR=$OUTDIR,INT_START=$INT_START,INT_STEP=$INT_STEP,N_INT=$N_INT,N_T=$N_T \
  array_job_no_imbalance.sh
