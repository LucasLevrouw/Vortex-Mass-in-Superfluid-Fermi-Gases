#!/bin/bash
# Run this to submit: bash launch_finT.sh

# --- Define ranges here ---
T_START=0.0;    T_STOP=0.95;  T_STEP=0.05
INT_START=-1.0; INT_STOP=1.0; INT_STEP=1.0
ZETA_START=0.0; ZETA_STOP=0.9; ZETA_STEP=0.1

# --- Auto-calculate counts ---
N_T=$(awk    "BEGIN { print int(($T_STOP    - $T_START)    / $T_STEP    + 0.5) + 1 }")
N_INT=$(awk  "BEGIN { print int(($INT_STOP  - $INT_START)  / $INT_STEP  + 0.5) + 1 }")
N_ZETA=$(awk "BEGIN { print int(($ZETA_STOP - $ZETA_START) / $ZETA_STEP + 0.5) + 1 }")
N_JOBS=$(( N_T * N_INT * N_ZETA ))

OUTDIR=out_finT

echo "Submitting $N_T x $N_INT x $N_ZETA = $N_JOBS jobs"

mkdir -p $OUTDIR

sbatch \
  --array=0-$(( N_JOBS - 1 )) \
  --output=/dev/null \
  --export=ALL,OUTDIR=$OUTDIR,T_START=$T_START,T_STEP=$T_STEP,N_T=$N_T,\
INT_START=$INT_START,INT_STEP=$INT_STEP,N_INT=$N_INT,\
ZETA_START=$ZETA_START,ZETA_STEP=$ZETA_STEP,N_ZETA=$N_ZETA \
  array_job_finT.sh