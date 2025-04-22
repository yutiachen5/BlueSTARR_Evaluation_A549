#!/bin/sh
#
#SBATCH --get-user-env
#SBATCH -A majoroslab
#SBATCH -J unbiased-pred
#SBATCH -o /datacommons/igvf-pm/A549/GR-AP1/simulated-seq/slurm/slurm-0.001/slurm-unbiased-train/outputs/%a.out
#SBATCH -e /datacommons/igvf-pm/A549/GR-AP1/simulated-seq/slurm/slurm-0.001/slurm-unbiased-train/outputs/%a.err
#SBATCH --exclusive
#SBATCH --gres=gpu:1
#SBATCH --array=1-5%5
#SBATCH -p majoroslab-gpu,igvf-gpu
#SBATCH --nice=100
#SBATCH --mem=102400
#SBATCH --cpus-per-task=1
#
/datacommons/igvf-pm/A549/GR-AP1/simulated-seq/slurm/slurm-0.001/slurm-unbiased-train/command${SLURM_ARRAY_TASK_ID}.sh
