# process raw data of A549 extra CGs
# process lines to have the format of >17698785 /coord=chr1:138900-139200
# CCTGTAGACGCTGACAGGAGGCAGGAGCTGGGCCTGGACAGGTCAACTTGAGGAGATTTT
# process counts to have the format
# DNA=5    RNA=4
# 1   2   3   4   5   6   7    8    9

# without normalization
# counts are rounded to the nearest integer

import pandas as pd
from FastaWriter import FastaWriter
from tqdm import tqdm
import sys
import gzip
import os
    
def main(input_dir,output_dir,labels):
    input_path = input_dir+'A549.'+label+'.combined.input_and_output.w300s50.log2FC.gt_200.sequence.final.txt.gz' # shreshold=200
    df_combined = pd.read_csv(input_path, sep = '\t', compression="gzip")

    output_path = output_dir+label+'/data/'
    os.makedirs(output_path, exist_ok=True)

    count_cols = ['input_rep1', 'input_rep2', 'input_rep3', 'input_rep4', 'input_rep5'] + \
                    ['A549_'+label+'_rep'+str(i) for i in range(2,6)]
    df_combined[count_cols] = df_combined[count_cols].round().astype(int)
    df_count = df_combined[count_cols]

    ls_chr = df_combined['chrom']
    ls_start = df_combined['start']
    ls_end = df_combined['end']
    ls_seq = df_combined['sequence']

    fastaWriter=FastaWriter()
    output_fasta_train = gzip.open(output_path+label+"-all-train.fasta.gz","wt") 
    output_fasta_val = gzip.open(output_path+label+"-all-val.fasta.gz","wt") 
    output_fasta_test = gzip.open(output_path+label+"-all-test.fasta.gz","wt") 
    output_count_train = gzip.open(output_path+label+"-all-train-counts.txt.gz","wt") 
    output_count_val = gzip.open(output_path+label+"-all-val-counts.txt.gz","wt") 
    output_count_test = gzip.open(output_path+label+"-all-test-counts.txt.gz","wt") 

    for i in tqdm(range(len(df_combined)), desc=label):
        defline = '>'+str(i)+' /coord='+str(ls_chr[i])+':'+str(ls_start[i])+'-'+str(ls_end[i])
        line = df_count.iloc[i]
        if ls_chr[i] == 'chr1':
            fastaWriter.addToFasta(defline,ls_seq[i],output_fasta_test) # chr1 for test
            print('\t'.join(map(str, line.values)),end="",file=output_count_train)
        elif ls_chr[i] == 'chr2':
            fastaWriter.addToFasta(defline,ls_seq[i],output_fasta_val) # chr2 for validation
            print('\t'.join(map(str, line.values)),end="",file=output_count_val)
        else:
            fastaWriter.addToFasta(defline,ls_seq[i],output_fasta_test)
            print('\t'.join(map(str, line.values)),end="",file=output_count_test)


    # df_count = df_combined[count_cols]

    # output_count = output_path+label+'-all-counts.txt.gz'
    # with open(output_count, 'w') as file:
    #     file.write("DNA=5\tRNA=4\n")  
    #     df_count.to_csv(file, sep='\t', index=False, header=None)  

if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <input-dir:/.../> <output-dir:/.../> <label:A,B,C>\n")
(input_dir,output_dir,labels)=sys.argv[1:]

labels = labels.split(',') if len(labels)>1 else labels.to_list()
for label in labels:
    main(input_dir,output_dir,label)

# command line:
# python /work/igvf-pm/A549/extra_GCs/process_combined.py /work/igvf-pm/alex_b/starrseq_A549_extra_GCs/ /work/igvf-pm/A549/extra_GCs/ AZD2906,AZD9567,CORT108297,CpdA,GW870086,Hydrocortisone,RU486,ZK216348