#!/bin/bash
#SBATCH --ntasks=1 --cpus-per-task=1
#SBATCH --time=20:00:00

module purge
module load calcua/2025a
module load Julia/1.12.4

i_int=$(( SLURM_ARRAY_TASK_ID / N_ZETA ))
i_zeta=$(( SLURM_ARRAY_TASK_ID % N_ZETA ))

int=$(awk  "BEGIN { printf \"%.2f\", $INT_START  + $i_int  * $INT_STEP  }")
zeta=$(awk "BEGIN { printf \"%.2f\", $ZETA_START + $i_zeta * $ZETA_STEP }")

exec > $OUTDIR/output_T=0_int=${int}_zeta=${zeta}.out 2>&1

julia --project=$SLURM_SUBMIT_DIR/.. $SLURM_SUBMIT_DIR/calculate-profile-one.jl 0 "$int" "$zeta"