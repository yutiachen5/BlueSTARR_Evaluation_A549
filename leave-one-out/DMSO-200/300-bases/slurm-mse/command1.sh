#!/bin/sh
hostname && nvidia-smi && python /datacommons/igvf-pm/K562/leave-one-out/BlueSTARR/BlueSTARR-multitask-aver.py /datacommons/igvf-pm/A549/leave-one-out/DMSO-200/config-mse/A549_chr1.config /datacommons/igvf-pm/A549/leave-one-out/DMSO-200/data /datacommons/igvf-pm/A549/leave-one-out/DMSO-200/slurm-mse/outputs/
