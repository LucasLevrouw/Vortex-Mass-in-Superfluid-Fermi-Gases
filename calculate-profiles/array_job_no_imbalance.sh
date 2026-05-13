#!/bin/bash
#SBATCH --ntasks=1 --cpus-per-task=1
#SBATCH --time=20:00:00

module purge
module load calcua/2025a
module load Julia/1.12.4

T_VALUES=(0.0 0.5 0.9)

i_T=$(( SLURM_ARRAY_TASK_ID / N_INT ))
i_int=$(( SLURM_ARRAY_TASK_ID % N_INT ))

T=${T_VALUES[$i_T]}
int=$(awk "BEGIN { printf \"%.2f\", $INT_START + $i_int * $INT_STEP }")
zeta="0.00"

exec > $OUTDIR/output_T=${T}_int=${int}_zeta=${zeta}.out 2>&1

julia --project=$SLURM_SUBMIT_DIR/.. $SLURM_SUBMIT_DIR/calculate-profile-one.jl "$T" "$int" "$zeta"
