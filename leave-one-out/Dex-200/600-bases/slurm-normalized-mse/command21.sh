#!/bin/sh
hostname && nvidia-smi
python /datacommons/igvf-pm/K562/leave-one-out/BlueSTARR/BlueSTARR-multitask-A549.py /datacommons/igvf-pm/A549/leave-one-out/Dex-200/config-mse/A549_chr21.config /datacommons/igvf-pm/A549/leave-one-out/Dex-200/600-bases/data-normalized /datacommons/igvf-pm/A549/leave-one-out/Dex-200/600-bases/slurm-normalized-mse/outputs/
