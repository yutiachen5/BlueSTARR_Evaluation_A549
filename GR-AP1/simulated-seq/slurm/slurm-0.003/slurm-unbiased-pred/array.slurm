#!/bin/sh
#
#SBATCH --get-user-env
#SBATCH -A majoroslab
#SBATCH -J pred
#SBATCH -o /datacommons/igvf-pm/A549/GR-AP1/simulated-seq/slurm/slurm-0.003/slurm-lognormal10-pred/outputs/%a.out
#SBATCH -e /datacommons/igvf-pm/A549/GR-AP1/simulated-seq/slurm/slurm-0.003/slurm-lognormal10-pred/outputs/%a.err
#SBATCH --exclusive
#SBATCH --gres=gpu:1
#SBATCH --array=1-5%5
#SBATCH -p majoroslab-gpu,igvf-gpu
#SBATCH --nice=100
#SBATCH --mem=102400
#SBATCH --cpus-per-task=1
#
/datacommons/igvf-pm/A549/GR-AP1/simulated-seq/slurm/slurm-0.003/slurm-lognormal10-pred/command${SLURM_ARRAY_TASK_ID}.sh
