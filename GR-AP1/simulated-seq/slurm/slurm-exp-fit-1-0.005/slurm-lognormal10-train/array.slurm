#!/bin/sh
#
#SBATCH --get-user-env
#SBATCH -A majoroslab
#SBATCH -J fit-1-train
#SBATCH -o /datacommons/igvf-pm/A549/GR-AP1/simulated-seq/slurm/slurm-exp-fit-1-0.005/slurm-lognormal10-train/outputs/%a.out
#SBATCH -e /datacommons/igvf-pm/A549/GR-AP1/simulated-seq/slurm/slurm-exp-fit-1-0.005/slurm-lognormal10-train/outputs/%a.err
#SBATCH --exclusive
#SBATCH --gres=gpu:1
#SBATCH --array=1-5%5
#SBATCH -p majoroslab-gpu,igvf-gpu
#SBATCH --nice=100
#SBATCH --mem=102400
#SBATCH --cpus-per-task=1
#
/datacommons/igvf-pm/A549/GR-AP1/simulated-seq/slurm/slurm-exp-fit-1-0.005/slurm-lognormal10-train/command${SLURM_ARRAY_TASK_ID}.sh
