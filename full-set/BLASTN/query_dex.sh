#!/bin/bash
#SBATCH -o /hpc/group/igvf/A549/full-set/BLASTN/query_dex.out
#SBATCH -e /hpc/group/igvf/A549/full-set/BLASTN/query_dex.err
#SBATCH -A majoroslab
#SBATCH -p majoroslab-gpu,igvf-gpu
#SBATCH --mem=10G

module load NCBI-BLAST

QUERY="/hpc/group/igvf/A549/processed-data/300-bases/Dex-200/Dex-200-all.fasta"
DATABASE="/hpc/group/igvf/A549/full-set/BLASTN/all_indexed_dex"
OUTPUT="/hpc/group/igvf/A549/full-set/BLASTN/all_aligned_300_dex.txt"
EVALUE="1e-5"
OUTFMT="6"


blastn -query $QUERY -db $DATABASE -out $OUTPUT -evalue $EVALUE -outfmt $OUTFMT -sum_stats true -perc_identity 90.0 -word_size 100
echo "Query completed. Results saved to "$OUTPUT