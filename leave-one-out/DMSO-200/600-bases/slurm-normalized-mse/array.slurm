#!/bin/sh
#
#SBATCH --get-user-env
#SBATCH -A majoroslab
#SBATCH -J A549
#SBATCH -o /datacommons/igvf-pm/A549/leave-one-out/DMSO-200/600-bases/slurm-normalized-mse/outputs/%a.out
#SBATCH -e /datacommons/igvf-pm/A549/leave-one-out/DMSO-200/600-bases/slurm-normalized-mse/outputs/%a.err
#SBATCH --exclusive
#SBATCH --gres=gpu:1
#SBATCH --array=1
#SBATCH -p majoroslab-gpu,allenlab-gpu,scavenger-gpu
#SBATCH --nice=100
#SBATCH --mem=102400
#SBATCH --cpus-per-task=1
#
/datacommons/igvf-pm/A549/leave-one-out/DMSO-200/600-bases/slurm-normalized-mse/command${SLURM_ARRAY_TASK_ID}.sh
