#!/bin/bash
#SBATCH -o /hpc/group/igvf/A549/full-set/BLASTN/build_dex.out
#SBATCH -e /hpc/group/igvf/A549/full-set/BLASTN/build_dex.err
#SBATCH -A majoroslab
#SBATCH -p majoroslab-gpu,igvf-gpu
#SBATCH --mem=10G

module load NCBI-BLAST

makeblastdb -in "/hpc/group/igvf/A549/processed-data/300-bases/Dex-200/Dex-200-all.fasta" -dbtype nucl -out "/hpc/group/igvf/A549/full-set/BLASTN/all_indexed_dex" -hash_index

echo "Custom database built."


