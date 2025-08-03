#!/bin/sh
hostname && nvidia-smi && python /datacommons/igvf-pm/K562/leave-one-out/BlueSTARR/BlueSTARR-multitask-aver.py /datacommons/igvf-pm/A549/leave-one-out/DMSO-200/config-custom/A549_chrX.config /datacommons/igvf-pm/A549/leave-one-out/DMSO-200/data /datacommons/igvf-pm/A549/leave-one-out/DMSO-200/slurm-custom/outputs/
