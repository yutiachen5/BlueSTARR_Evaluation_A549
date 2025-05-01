#!/bin/sh
#
#SBATCH --get-user-env
#SBATCH -A majoroslab
#SBATCH -J cliff-0.5
#SBATCH -o /hpc/group/igvf/A549/GR-AP1/simulated-seq/slurm/cliff/no_enrichment_0.5/pred/outputs/%a.out
#SBATCH -e /hpc/group/igvf/A549/GR-AP1/simulated-seq/slurm/cliff/no_enrichment_0.5/pred/outputs/%a.err
#SBATCH --exclusive
#SBATCH --gres=gpu:1
#SBATCH --array=1-5%5
#SBATCH -p gpu-common,scavenger-gpu,biostat-gpu,majoroslab-gpu,igvf-gpu
#SBATCH --nice=100
#SBATCH --mem=102400
#SBATCH --cpus-per-task=1
#
/hpc/group/igvf/A549/GR-AP1/simulated-seq/slurm/cliff/no_enrichment_0.5/pred/command${SLURM_ARRAY_TASK_ID}.sh
