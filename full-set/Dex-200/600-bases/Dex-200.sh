#!/bin/sh
#
#SBATCH --get-user-env
#SBATCH -J dex-full
#SBATCH -A majoroslab
#SBATCH -o /datacommons/igvf-pm/A549/full-set/Dex-200/600-bases/outputs/Dex-200-1010.out
#SBATCH -e /datacommons/igvf-pm/A549/full-set/Dex-200/600-bases/outputs/Dex-200-1010.err
#SBATCH --exclusive
#SBATCH --gres=gpu:1
#SBATCH -p scavenger-gpu,majoroslab-gpu,igvf-gpu
#SBATCH --nice=100
#SBATCH --mem=102400
#SBATCH --cpus-per-task=1
#

hostname 
nvidia-smi 
python /datacommons/igvf-pm/A549/full-set/BlueSTARR/BlueSTARR-multitask.py /datacommons/igvf-pm/A549/full-set/BlueSTARR/A549-600-Dex-1010.config /datacommons/igvf-pm/A549/full-set/Dex-200/600-bases/data-normalized /datacommons/igvf-pm/A549/full-set/Dex-200/600-bases/saved_models/Dex-200-1010
