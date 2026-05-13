#!/bin/bash
# Run this to submit: bash launch_T0.sh

# --- Define ranges here ---
INT_START=-1.0;  INT_STOP=1.0;   INT_STEP=0.05
ZETA_START=0.0;  ZETA_STOP=0.95; ZETA_STEP=0.05

# --- Auto-calculate counts ---
N_INT=$(awk  "BEGIN { print int(($INT_STOP  - $INT_START)  / $INT_STEP + 0.5) + 1 }")
N_ZETA=$(awk "BEGIN { print int(($ZETA_STOP - $ZETA_START) / $ZETA_STEP + 0.5) + 1 }")
N_JOBS=$(( N_INT * N_ZETA ))

OUTDIR=out_T0

echo "Submitting $N_INT x $N_ZETA = $N_JOBS jobs"

mkdir -p $OUTDIR

sbatch \
  --array=0-$(( N_JOBS - 1 )) \
  --output=/dev/null \
  --export=ALL,OUTDIR=$OUTDIR,INT_START=$INT_START,INT_STEP=$INT_STEP,N_INT=$N_INT,\
ZETA_START=$ZETA_START,ZETA_STEP=$ZETA_STEP,N_ZETA=$N_ZETA \
  array_job_T0.sh